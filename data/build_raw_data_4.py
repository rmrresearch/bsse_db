"""
build_raw_data_4.py
--------------------
Parses all NWChem .out files for the H2O tetramer from 4.tar.gz and stores
the extracted energies in a nested dictionary mirroring the directory structure.

Path structure inside the tarball:
    4 / CCSD_T / aug-cc-pvdz / REAL_BASIS.out
    [0]   [1]       [2]          stem

Since the tetramer has a single geometry (no distance scan), config is
hardcoded as '0'.

Dictionary structure
--------------------
    raw_data['H2O']['4']['0']['SCF']['aug-cc-pvdz']['A_AB']    = E_scf
    raw_data['H2O']['4']['0']['MP2']['aug-cc-pvdz']['A_AB']    = mp2_corr
    raw_data['H2O']['4']['0']['CCSD_T']['aug-cc-pvdz']['A_AB'] = ccsdt_comp

    parameters['H2O']['4']['0'] = {
        'config_label': '0',
        'params': {}
    }

Note on stored values
---------------------
    'SCF'    : total SCF energy (Hartree)
    'MP2'    : MP2 correlation energy  = E_mp2_total   - E_scf
    'CCSD_T' : CCSD(T) component       = E_ccsdt_total - E_mp2_total

To recover total energies:
    E_scf_total   = raw_data[...]['SCF'][basis][key]
    E_mp2_total   = E_scf_total + raw_data[...]['MP2'][basis][key]
    E_ccsdt_total = E_mp2_total + raw_data[...]['CCSD_T'][basis][key]

Usage
-----
    python3 build_raw_data_4.py
    # produces data/raw_data_H2O_4.pickle

    from build_raw_data_4 import build_raw_data
    raw_data, raw_data_noconv, parameters = build_raw_data('H2O')
"""

import re
import pickle
from pathlib import Path, PurePosixPath

from tar_output_discovery_4 import iter_out_texts


# ---------------------------------------------------------------------------
# Energy extraction
# ---------------------------------------------------------------------------

# NWChem output patterns
_SCF_PAT   = re.compile(r'Total SCF energy\s*[=:]\s*([-\d.]+)')
_RHF_PAT   = re.compile(r'Total RHF energy\s*[=:]\s*([-\d.]+)')
_MP2_PAT   = re.compile(r'Total MP2 energy\s*[=:]\s*([-\d.]+)')
_CCSD_PAT  = re.compile(r'CCSD\(T\) total energy\s*[=:]\s*([-\d.]+)')
_CCSD_PAT2 = re.compile(r'Total CCSD\(T\) energy\s*[=:]\s*([-\d.]+)')


def parse_energies(text):
    """
    Extract SCF, MP2, and CCSD(T) total energies from a NWChem output file.

    Returns
    -------
    (E_scf, E_mp2, E_ccsdt) : floats or None if not found
    """
    E_scf   = None
    E_mp2   = None
    E_ccsdt = None

    for line in text.splitlines():
        if E_scf is None:
            m = _SCF_PAT.search(line) or _RHF_PAT.search(line)
            if m:
                E_scf = float(m.group(1))
        if E_mp2 is None:
            m = _MP2_PAT.search(line)
            if m:
                E_mp2 = float(m.group(1))
        if E_ccsdt is None:
            m = _CCSD_PAT.search(line) or _CCSD_PAT2.search(line)
            if m:
                E_ccsdt = float(m.group(1))

    return E_scf, E_mp2, E_ccsdt


# ---------------------------------------------------------------------------
# Main builder
# ---------------------------------------------------------------------------

def build_raw_data(molecule, data_dir=None):
    """
    Parse all tetramer .out files and return nested energy dictionaries.

    Parameters
    ----------
    molecule : str
        Molecule name, e.g. 'H2O'
    data_dir : Path or None
        Path to the molecule directory (e.g. data/H2O/).
        Defaults to <script_dir>/<molecule>/

    Returns
    -------
    raw_data       : nested dict  — parsed energies (correlation/component form)
    raw_data_noconv: nested dict  — files that failed to converge
    parameters     : nested dict  — config metadata
    """
    if data_dir is None:
        here     = Path(__file__).resolve().parent
        data_dir = here / molecule

    raw_data        = {}
    raw_data_noconv = {}
    parameters      = {}

    cluster = '4'
    config  = '0'                     # single geometry — no distance scan

    for tar_path, internal_path, text in iter_out_texts(data_dir):

        # Only process tetramer files
        if internal_path.parts[0] != cluster:
            continue

        # Unpack path:  4 / CCSD_T / aug-cc-pvdz / REAL_BASIS.out
        parts      = internal_path.parts
        method_dir = parts[1]           # 'CCSD_T'
        basis      = parts[2]           # 'aug-cc-pvdz'
        file_key   = PurePosixPath(parts[3]).stem   # e.g. 'A_AB'

        # Extract energies
        E_scf, E_mp2, E_ccsdt = parse_energies(text)

        # Track non-converged files
        if E_scf is None or E_mp2 is None or E_ccsdt is None:
            bucket = raw_data_noconv
            target = _nested(bucket, molecule, cluster, config,
                             basis, file_key)
            target['SCF']    = E_scf
            target['MP2']    = E_mp2
            target['CCSD_T'] = E_ccsdt
            continue

        # Compute correlation / component energies
        mp2_corr      = E_mp2   - E_scf
        ccsdt_comp    = E_ccsdt - E_mp2

        # Store in raw_data
        for method, value in [('SCF',    E_scf),
                               ('MP2',    mp2_corr),
                               ('CCSD_T', ccsdt_comp)]:
            d = raw_data
            for key in [molecule, cluster, config, method, basis]:
                d = d.setdefault(key, {})
            d[file_key] = value

    # Record config metadata
    parameters \
        .setdefault(molecule, {}) \
        .setdefault(cluster, {})[config] = {
            'config_label': config,
            'params': {}
        }

    return raw_data, raw_data_noconv, parameters


def _nested(d, *keys):
    """Navigate/create nested dict levels, return the innermost dict."""
    for key in keys:
        d = d.setdefault(key, {})
    return d


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

if __name__ == '__main__':
    here     = Path(__file__).resolve().parent
    data_dir = here / 'H2O'
    molecule = 'H2O'

    print('Parsing tetramer output files...')
    raw_data, raw_data_noconv, parameters = build_raw_data(
        molecule, data_dir=data_dir
    )

    # --- Convergence report ---
    n_keys = len(raw_data.get(molecule, {})
                         .get('4', {})
                         .get('0', {})
                         .get('SCF', {})
                         .get('aug-cc-pvdz', {}))
    print(f'Converged file keys : {n_keys}  (expected 65)')

    if raw_data_noconv:
        print(f'WARNING — {len(raw_data_noconv)} non-converged file(s):')
        for k in raw_data_noconv:
            print(f'  {k}')
    else:
        print('All files converged successfully. ✅')

    # --- Sample checks ---
    basis = 'aug-cc-pvdz'
    print('\nSample values:')
    for key in ['A_A', 'ABCD_ABCD', 'A_ABCD']:
        scf = raw_data[molecule]['4']['0']['SCF'][basis].get(key)
        print(f"  SCF  [{key}] = {scf}")

    print(f'\nConfig metadata : {parameters[molecule]["4"]}')

    # --- Pickle ---
    bundle = {
        'raw_data'  : raw_data,
        'parameters': parameters,
    }
    out_path = here / f'raw_data_{molecule}_4.pickle'
    with open(out_path, 'wb') as f:
        pickle.dump(bundle, f)
    print(f'\nPickled to: {out_path}')