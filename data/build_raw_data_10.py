"""
build_raw_data_10.py
--------------------
Parses all NWChem .out files for the HF decamer directly from disk
(no tar navigation) and stores extracted energies in a nested dictionary
mirroring the directory structure.

Path structure (flat directory, 4135 .out files):
    HF/10/CCSD_T/aug-cc-pvdz/REAL_BASIS.out

Single geometry (config '0'), single basis (aug-cc-pvdz).
All 10 monomers A-J are geometrically unique -- no injection needed.
Calculations are full CCSD(T)/aug-cc-pvdz, truncated at 4-body basis sets.

Dictionary structure
--------------------
    raw_data['HF']['10']['0']['SCF']['aug-cc-pvdz']['A_ABCD']    = E_scf
    raw_data['HF']['10']['0']['MP2']['aug-cc-pvdz']['A_ABCD']    = mp2_corr
    raw_data['HF']['10']['0']['CCSD_T']['aug-cc-pvdz']['A_ABCD'] = ccsdt_comp

    parameters['HF']['10']['0'] = {
        'config_label': '0',
        'params': {}
    }

Stored values (incremental, not cumulative)
-------------------------------------------
    'SCF'    : total SCF energy (Hartree)
    'MP2'    : MP2 correlation energy  = E_mp2_total - E_scf
    'CCSD_T' : CCSD(T) component       = E_ccsdt_total - E_mp2_total

Expected key counts by body order of real fragment
---------------------------------------------------
    1-body real  : 10 + 90 + 360 + 840  = 1300  keys
    2-body real  : 45 + 360 + 1260      = 1665  keys
    3-body real  : 120 + 840            =  960  keys
    4-body real  : 210                  =  210  keys
    Total                               = 4135  keys

Usage
-----
    python3 build_raw_data_10.py
    # produces  data/raw_data_HF_10.pickle
"""

import re
import pickle
from pathlib import Path


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

MOLECULE       = 'HF'
CLUSTER        = '10'
CONFIG         = '0'
BASIS          = 'aug-cc-pvdz'
METHOD_DIR     = 'CCSD_T'
EXPECTED_FILES = 4135

# Expected key counts by body order of real fragment (for verification)
EXPECTED_BY_BODY = {1: 1300, 2: 1665, 3: 960, 4: 210}


# ---------------------------------------------------------------------------
# Energy extraction  (same patterns as build_raw_data_4.py)
# ---------------------------------------------------------------------------

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
    (E_scf, E_mp2, E_ccsdt) : floats, or None if not found
    """
    E_scf = E_mp2 = E_ccsdt = None

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

def build_raw_data(data_dir=None):
    """
    Parse all HF decamer .out files and return nested energy dictionaries.

    Parameters
    ----------
    data_dir : Path or None
        Path to the HF molecule directory (e.g. data/HF/).
        Defaults to <script_dir>/HF/

    Returns
    -------
    raw_data        : nested dict  -- parsed energies (incremental form)
    raw_data_noconv : nested dict  -- files that failed to converge
    parameters      : nested dict  -- config metadata
    n_parsed        : int          -- number of successfully parsed files
    """
    here = Path(__file__).resolve().parent

    if data_dir is None:
        data_dir = here / MOLECULE

    out_dir = data_dir / CLUSTER / METHOD_DIR / BASIS

    if not out_dir.exists():
        raise FileNotFoundError(
            f'Output directory not found: {out_dir}\n'
            f'Expected: data/HF/10/CCSD_T/aug-cc-pvdz/'
        )

    raw_data        = {}
    raw_data_noconv = {}
    parameters      = {}
    n_parsed        = 0

    out_files = sorted(out_dir.glob('*.out'))

    for out_file in out_files:
        file_key = out_file.stem          # e.g. 'A_ABCD'
        text     = out_file.read_text()

        E_scf, E_mp2, E_ccsdt = parse_energies(text)

        # Track non-converged or incomplete files
        if E_scf is None or E_mp2 is None or E_ccsdt is None:
            (raw_data_noconv
             .setdefault(MOLECULE, {})
             .setdefault(CLUSTER, {})
             .setdefault(CONFIG, {})
             .setdefault(BASIS, {}))[file_key] = {
                'SCF'   : E_scf,
                'MP2'   : E_mp2,
                'CCSD_T': E_ccsdt,
            }
            continue

        # Incremental components
        mp2_corr   = E_mp2   - E_scf
        ccsdt_comp = E_ccsdt - E_mp2

        # Store all three levels
        for method, value in [('SCF',    E_scf),
                               ('MP2',    mp2_corr),
                               ('CCSD_T', ccsdt_comp)]:
            (raw_data
             .setdefault(MOLECULE, {})
             .setdefault(CLUSTER, {})
             .setdefault(CONFIG, {})
             .setdefault(method, {})
             .setdefault(BASIS, {}))[file_key] = value

        n_parsed += 1

    # Record config metadata
    (parameters
     .setdefault(MOLECULE, {})
     .setdefault(CLUSTER, {}))[CONFIG] = {
        'config_label': CONFIG,
        'params': {}
    }

    return raw_data, raw_data_noconv, parameters, n_parsed


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

if __name__ == '__main__':
    here = Path(__file__).resolve().parent

    print('Parsing HF decamer output files...')
    print(f'  Source: HF/10/CCSD_T/aug-cc-pvdz/\n')

    raw_data, raw_data_noconv, parameters, n_parsed = build_raw_data()

    scf_keys = (raw_data
                .get(MOLECULE, {})
                .get(CLUSTER, {})
                .get(CONFIG, {})
                .get('SCF', {})
                .get(BASIS, {}))

    # --- Convergence report ---
    n_keys = len(scf_keys)
    status = '✅' if n_keys == EXPECTED_FILES else f'⚠️  expected {EXPECTED_FILES}'
    print(f'Converged file keys : {n_keys}  {status}')

    if raw_data_noconv:
        flat = (raw_data_noconv
                .get(MOLECULE, {})
                .get(CLUSTER, {})
                .get(CONFIG, {})
                .get(BASIS, {}))
        print(f'WARNING -- {len(flat)} non-converged file(s):')
        for k in sorted(flat):
            d = flat[k]
            print(f'  {k}  SCF={d["SCF"]}  MP2={d["MP2"]}  CCSD_T={d["CCSD_T"]}')
    else:
        print('All files converged successfully. ✅')

    # --- Key counts by body order of real fragment ---
    print('\nKey counts by body order of real fragment:')
    all_ok = True
    for n_body, expected in EXPECTED_BY_BODY.items():
        count = sum(
            1 for k in scf_keys
            if len(k.split('_')[0]) == n_body
        )
        ok = count == expected
        flag = '✅' if ok else f'⚠️  expected {expected}'
        print(f'  {n_body}-body real: {count:5d}  {flag}')
        if not ok:
            all_ok = False

    if all_ok:
        print('  All body-order counts correct. ✅')

    # --- Sample SCF values ---
    print('\nSample SCF values (Ha):')
    sample_keys = [
        'A_A', 'B_B', 'J_J',         # monomers in own basis
        'AB_AB', 'IJ_IJ',             # dimers in own basis
        'A_AB', 'B_AB',               # monomers in dimer basis
        'ABC_ABC',                     # trimer in own basis
        'A_ABC', 'AB_ABC',            # sub-clusters in trimer basis
        'ABCD_ABCD',                  # quadruple in own basis
        'A_ABCD', 'AB_ABCD',          # sub-clusters in quad basis
        'ABC_ABCD',                   # trimer in quad basis
    ]
    for key in sample_keys:
        val = scf_keys.get(key)
        if val is not None:
            print(f'  [{key:15s}] = {val:+.10f} Ha')
        else:
            print(f'  [{key:15s}] = None  ⚠️')

    # --- Incremental component check on one monomer ---
    print('\nIncremental component check for A_A:')
    for method in ['SCF', 'MP2', 'CCSD_T']:
        val = (raw_data
               .get(MOLECULE, {})
               .get(CLUSTER, {})
               .get(CONFIG, {})
               .get(method, {})
               .get(BASIS, {})
               .get('A_A'))
        label = {
            'SCF'   : 'E_scf              ',
            'MP2'   : 'E_mp2 - E_scf      ',
            'CCSD_T': 'E_ccsdt - E_mp2    ',
        }[method]
        print(f'  {method:8s} ({label}) = {val:+.10f} Ha')

    print(f'\nConfig metadata: {parameters[MOLECULE][CLUSTER]}')

    # --- Pickle ---
    out_path = here / f'raw_data_{MOLECULE}_{CLUSTER}.pickle'
    with open(out_path, 'wb') as f:
        pickle.dump({'raw_data': raw_data, 'parameters': parameters}, f)
    print(f'\nPickled to: {out_path}')