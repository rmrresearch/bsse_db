"""
build_raw_data_3.py  (generalized)
====================================
Parser for trimer (cluster size 3) NWChem output files.
Handles H2O, Ne, and HF trimers.

Reads the canonical tarball  data/<molecule>/3.tar.gz  without extracting it
to disk. For each .out file it extracts three total energies:
    - Total SCF energy
    - Total MP2 energy
    - Total CCSD(T) energy

From the total energies it computes:
    - MP2 correlation energy        : E_mp2   - E_scf
    - CCSD(T) correlation component : E_ccsdt - E_mp2

The parsed data is stored in a nested dictionary (raw_data) mirroring the
directory structure:

    raw_data[molecule]['3'][config_key][level_of_theory][basis][file_key]

For example:
    raw_data['Ne']['3']['0']['CCSD_T']['aug-cc-pvdz']['A_AB'] = energy

A separate parameters dictionary stores the geometry information:

    parameters[molecule]['3'][config_key] = {
        'config_label': '0_1.5',
        'params': {'distance': 1.5}
    }

Configurations missing any of the three energies are excluded from raw_data
and recorded in raw_data_noconv for inspection.

The final raw_data and parameters are bundled and pickled to
data/raw_data_{molecule}_3.pickle via pickle_data.py.

Notes
-----
- All 19 file keys are present for all molecules (confirmed for Ne, HF, H2O).
  B_B and C_C are stored as separate files — no aliasing needed.
- This script filters by parts[0] == '3' so only trimer files are processed
  even if 2.tar.gz is also present.
- This script lives in data/ alongside all other parsers.

Usage
-----
Run from anywhere; all paths are resolved relative to the script location:

    python build_raw_data_3.py --mol Ne
    python build_raw_data_3.py --mol HF
    python build_raw_data_3.py --mol H2O

Or interactively from the console (run from data/):

    from build_raw_data_3 import build_raw_data
    raw_data, raw_data_noconv, parameters = build_raw_data('Ne')
"""

import re
import argparse
from pathlib import Path, PurePosixPath
from tar_output_discovery_3 import iter_out_texts

# --- Regex patterns ---
SCF_PATTERN   = re.compile(r"Total\s+SCF\s+energy\s*[:=]\s*([\-0-9.+Ee]+)")
MP2_PATTERN   = re.compile(r"Total\s+MP2\s+energy\s*[:=]\s*([\-0-9.+Ee]+)")
CCSDT_PATTERN = re.compile(r"Total\s+CCSD\(T\)\s+energy\s*[:=]\s*([\-0-9.+Ee]+)")


def extract_last(pattern, text):
    """Return the last float matched by pattern in text, or None."""
    value = None
    for line in text.splitlines():
        m = pattern.search(line)
        if m:
            try:
                value = float(m.group(1))
            except ValueError:
                pass
    return value


def ensure_path(d, keys):
    """Ensure nested dict structure exists and return the deepest dict."""
    cur = d
    for k in keys:
        cur = cur.setdefault(k, {})
    return cur


def missing_energies(E_scf, E_mp2, E_ccsdt):
    """Return list of missing energy labels for a given output file."""
    missing = []
    if E_scf   is None: missing.append("SCF")
    if E_mp2   is None: missing.append("MP2")
    if E_ccsdt is None: missing.append("CCSD(T)")
    return missing


def build_raw_data(molecule, data_dir=None, verbose=True):
    """
    Build raw_data and parameters dictionaries for the trimer (cluster size 3).

    Works for H2O, Ne, and HF trimers. All 19 file keys are expected to be
    present as separate .out files for all molecules.

    raw_data is keyed by:
      raw_data[molecule]['3'][config_key][method_key][basis][file_key] = energy

    parameters is keyed by:
      parameters[molecule]['3'][config_key] = {
          'config_label': '0_1.5',
          'params': {'distance': 1.5}
      }

    Parameters
    ----------
    molecule : str
        Name of the molecule — 'H2O', 'Ne', or 'HF'.
    data_dir : Path, optional
        Path to the molecule data directory (e.g., .../data/Ne/).
        If not provided, derived from the script's location.
    verbose : bool
        If True, print convergence summary.
    """
    # If no data_dir provided, derive it from the script's location:
    # __file__ = data/build_raw_data_3.py
    # .parent / molecule = data/Ne/  (or data/HF/, data/H2O/)
    if data_dir is None:
        data_dir = Path(__file__).resolve().parent / molecule

    raw_data        = {}
    raw_data_noconv = {}
    parameters      = {}

    bucket = {}
    bad    = set()
    config_label_by_cfg = {}  # config_key -> config_label

    for tar_path, internal_path, text in iter_out_texts(data_dir):
        internal_path = PurePosixPath(str(internal_path))

        # Only process trimer files (parts[0] == '3')
        if len(internal_path.parts) < 5 or internal_path.parts[0] != '3':
            continue

        # Unpack path structure:
        # 3 / 0_1.5 / CCSD_T / aug-cc-pvdz / A_AB.out
        cluster      = internal_path.parts[0]          # '3'
        config_label = internal_path.parts[1]          # '0_1.5'
        method_dir   = internal_path.parts[2]          # 'CCSD_T'
        basis        = internal_path.parts[3]          # 'aug-cc-pvdz'
        file_key     = internal_path.stem              # 'A_AB'

        config_key = config_label.split("_")[0]        # '0'
        config_label_by_cfg[config_key] = config_label

        basis_folder_id = (cluster, config_key, method_dir, basis)

        # Extract energies
        E_scf   = extract_last(SCF_PATTERN,   text)
        E_mp2   = extract_last(MP2_PATTERN,   text)
        E_ccsdt = extract_last(CCSDT_PATTERN, text)

        if (E_scf is None) or (E_mp2 is None) or (E_ccsdt is None):
            bad.add(basis_folder_id)

            miss    = missing_energies(E_scf, E_mp2, E_ccsdt)
            cfg_bad = ensure_path(raw_data_noconv, (molecule, cluster, config_key))
            cfg_bad.setdefault("failures", []).append({
                "missing":       miss,
                "method_dir":    method_dir,
                "basis":         basis,
                "file_key":      file_key,
                "tar_path":      str(tar_path),
                "internal_path": str(internal_path),
            })
            continue

        mp2_corr        = E_mp2   - E_scf
        ccsdt_component = E_ccsdt - E_mp2

        bucket.setdefault(basis_folder_id, []).append(
            (file_key, E_scf, mp2_corr, ccsdt_component)
        )

    # --- Commit only good folders ---
    converged_cfgs = set()

    for (cluster, config_key, method_dir, basis), records in bucket.items():
        if (cluster, config_key, method_dir, basis) in bad:
            continue

        converged_cfgs.add(config_key)

        for file_key, E_scf, mp2_corr, ccsdt_component in records:
            cfg_node = ensure_path(raw_data, (molecule, cluster, config_key))
            ensure_path(cfg_node, ("SCF",    basis))[file_key] = E_scf
            ensure_path(cfg_node, ("MP2",    basis))[file_key] = mp2_corr
            ensure_path(cfg_node, ("CCSD_T", basis))[file_key] = ccsdt_component

    # --- Build parameters for converged configs only ---
    for config_key in converged_cfgs:
        config_label = config_label_by_cfg.get(config_key)
        if config_label is None:
            continue
        # Config label format is n_R (e.g. '0_1.5', '5_2.75')
        # Distance is always the second part after splitting on '_'
        distance = float(config_label.split("_")[1])
        ensure_path(parameters, (molecule, '3'))[config_key] = {
            "config_label": config_label,
            "params":       {"distance": distance},
        }

    if verbose:
        n_converged = len(converged_cfgs)
        print(f"\n{'='*50}")
        print(f"  Molecule : {molecule}   Cluster : 3")
        print(f"  Converged configs : {n_converged}")
        if raw_data_noconv:
            print("\n⚠️  WARNING: Some configurations did NOT converge.")
            print("    Please inspect `raw_data_noconv` for details.\n")
            for mol, mol_data in raw_data_noconv.items():
                for cluster, cluster_data in mol_data.items():
                    bad_cfgs = sorted(cluster_data.keys(), key=int)
                    print(f"    {mol} | cluster {cluster} | failed configs: {bad_cfgs}")
        else:
            print("  ✅ All configurations converged successfully.")
        print(f"{'='*50}\n")

    return raw_data, raw_data_noconv, parameters


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Parse trimer NWChem output files for a given molecule."
    )
    parser.add_argument(
        "--mol",
        type=str,
        required=True,
        choices=["H2O", "Ne", "HF"],
        help="Molecule name: H2O, Ne, or HF",
    )
    args = parser.parse_args()
    molecule = args.mol

    here     = Path(__file__).resolve().parent   # data/
    data_dir = here / molecule                   # data/Ne/, data/HF/, data/H2O/
    root_dir = here                              # data/

    raw_data, raw_data_noconv, parameters = build_raw_data(
        molecule = molecule,
        data_dir = data_dir,
    )

    from pickle_data import save_raw_bundle
    save_raw_bundle(
        molecule   = molecule,
        raw_data   = raw_data,
        parameters = parameters,
        outdir     = root_dir,
        cluster    = "3",
    )

    print(f"Pickled raw_data to {root_dir / f'raw_data_{molecule}_3.pickle'}")