# build_raw_data.py
"""
build_raw_data.py
=================
Parser for NWChem CCSD(T) output files produced in the BSSE database project
(rmrresearch/bsse_db, branch fcuantico, Issue #46).

Purpose
-------
This script implements Task 3 of the PI task list: dynamically discover all
NWChem .out files stored in compressed tarballs, extract the relevant energies,
and store them in a structured nested dictionary that mirrors the directory
hierarchy. The result is pickled to disk for downstream use by the many-body
interaction and BSSE computation scripts (Tasks 4 and 5).

Supported molecules
-------------------
H2O  — water
         dimer (2): translational dataset, aug-cc-pVDZ, 268 configs
         trimer (3): translational + rotational datasets
         tetramer (4): translational + rotational datasets
Ne   — neon
         dimer (2): translational dataset, aug-cc-pVDZ, 268 configs
                    radial dataset,        aug-cc-pVTZ,   5 configs
HF   — hydrogen fluoride
         dimer (2): translational dataset, aug-cc-pVDZ, 268 configs
                    radial dataset,        aug-cc-pVTZ,   3 configs
                    rotational dataset,    aug-cc-pVDZ, 152 configs

Directory structure expected inside each tarball
-------------------------------------------------
Each tarball is expected to contain paths of the form:

    <cluster>/<motion>/<config_label>/<METHOD>/<basis>/<file_key>.out

where:
    cluster      — cluster size integer as string, e.g. "2" for dimer
    motion       — one of: translational, radial, rotational (or rotational*)
    config_label — configuration folder, format n_<params>, e.g. "0_0.25_0.25_1.5"
                   n  : integer configuration index (lexicographic order)
                   params : x_y_z for trnl/radial; alpha_beta_gamma_R for rot
    METHOD       — level of theory directory, e.g. "CCSD_T"
    basis        — basis set directory, e.g. "aug-cc-pvdz" or "aug-cc-pvtz"
    file_key     — one of: A_A, AB_AB, A_AB, B_AB
                   following the <real>_<basis_set_monomers>.out convention
                   defined in Issue #46.

File key conventions (dimer)
-----------------------------
    A_A    : monomer A energy computed in monomer A basis
    AB_AB  : dimer AB energy computed in dimer AB basis
    A_AB   : monomer A energy computed in dimer AB basis (ghost functions on B)
    B_AB   : monomer B energy computed in dimer AB basis (ghost functions on A)

Note: for Ne and HF dimers, monomers A and B have identical geometries, so
no separate B_B calculation exists. The BSSE computation script uses A_A for
both monomers.

Energies extracted per .out file
---------------------------------
For each .out file the script extracts (using the last occurrence in the file):
    E_scf   : total SCF energy          (Hartree)
    E_mp2   : total MP2 energy          (Hartree)
    E_ccsdt : total CCSD(T) energy      (Hartree)

From these it computes and stores:
    SCF energy             = E_scf
    MP2 correlation energy = E_mp2   - E_scf
    CCSD(T) component      = E_ccsdt - E_mp2

Output data structures
----------------------
raw_data[molecule][cluster][motion][config_id][method_key][basis][file_key]
    Innermost value is a float energy (Ha).
    method_key is one of: "SCF", "MP2", "CCSD_T"

raw_data_noconv[molecule][cluster][motion][config_id]["failures"]
    List of dicts recording which energies were missing, which file and
    tarball they came from. Configurations with any missing energy are
    excluded from raw_data entirely.

parameters[molecule][cluster][motion][config_id]
    Dict with keys:
        "config_label" : full config folder name, e.g. "0_0.25_0.25_1.5"
        "params"       : dict of float geometric parameters, e.g.
                         {"x": 0.25, "y": 0.25, "z": 1.5}  for trnl/radial
                         {"alpha": 0.0, "beta": 90.0, "gamma": 0.0, "R": 2.75} for rot
    Only populated for configurations that appear in raw_data (converged).

Output pickle
-------------
One pickle file per molecule: raw_data_<MOL>.pickle
Each pickle contains a single dict with three keys:
    "raw_data"        -> raw_data dict
    "raw_data_noconv" -> raw_data_noconv dict
    "parameters"      -> parameters dict

Usage
-----
Run from bsse_db/data/ with the tarballs present, e.g.:

    # single molecule
    python3 build_raw_data.py --mol Ne

    # two molecules
    python3 build_raw_data.py --mol Ne HF

    # all three
    python3 build_raw_data.py --mol H2O Ne HF

Dependencies
------------
    tar_output_discovery.py  — must be in the same directory; provides
                               iter_out_texts(molecule) which walks all
                               *.tar.gz files under the molecule directory
                               and yields (tar_path, internal_path, text)
                               for every .out file found.

Notes
-----
- Python reads tarballs directly without extracting to disk.
- The parser is designed to be re-run as new tarballs are added; it will
  simply overwrite the pickle with updated data.
- The "last occurrence" strategy for energy extraction handles NWChem output
  files where energies are printed multiple times during iterative convergence;
  the final value is always the converged result.
"""

import re
import pickle
import argparse
from pathlib import PurePosixPath
from tar_output_discovery import iter_out_texts

# ---------------------------------------------------------------------------
# Regex patterns
# ---------------------------------------------------------------------------
SCF_PATTERN   = re.compile(r"Total\s+SCF\s+energy\s*[:=]\s*([\-0-9.+Ee]+)")
MP2_PATTERN   = re.compile(r"Total\s+MP2\s+energy\s*[:=]\s*([\-0-9.+Ee]+)")
CCSDT_PATTERN = re.compile(r"Total\s+CCSD\(T\)\s+energy\s*[:=]\s*([\-0-9.+Ee]+)")

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def extract_last(pattern, text):
    """
    Return the last float matched by *pattern* in *text*, or None.

    Parameters
    ----------
    pattern : re.Pattern
        Compiled regular expression with one capturing group for the numeric
        value, e.g. SCF_PATTERN, MP2_PATTERN, or CCSDT_PATTERN.
    text : str
        Full text content of a NWChem .out file.

    Returns
    -------
    float or None
        The last numeric value matched in the file, or None if no match was
        found. Using the last occurrence handles NWChem output files where
        energies are printed at each SCF/correlation iteration; the final
        printed value is always the converged result.
    """
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
    """
    Ensure a chain of nested dicts exists and return the deepest one.

    Parameters
    ----------
    d : dict
        The root dictionary to traverse or extend.
    keys : iterable
        Sequence of keys defining the path, e.g. (molecule, cluster, motion).

    Returns
    -------
    dict
        The innermost dict at the end of the key chain. Creates intermediate
        dicts with setdefault if they do not already exist.

    Example
    -------
    >>> d = {}
    >>> node = ensure_path(d, ("Ne", "2", "trnl"))
    >>> node["0"] = {}
    >>> d
    {'Ne': {'2': {'trnl': {'0': {}}}}}
    """
    cur = d
    for k in keys:
        cur = cur.setdefault(k, {})
    return cur


def motion_label(component):
    """
    Map the second path component to a short motion label.

      translational  -> trnl
      rotational*    -> rot
      radial         -> radial   (kept as-is)
      anything else  -> kept as-is
    """
    c = component.lower()
    if c.startswith("trans"):
        return "trnl"
    if c.startswith("rot"):
        return "rot"
    return component          # covers "radial" and any future motion types


def parse_params(motion, config_label):
    """
    Parse geometric parameters from a config_label string.

    Supported motions and label formats
    ------------------------------------
    trnl   : n_x_y_z              e.g. "0_0.25_0.25_1.5"
    radial : n_x_y_z              e.g. "0_0.0_0.0_2.5"   (same format as trnl)
    rot    : n_alpha_beta_gamma_R  e.g. "3_0.0_90.0_0.0_2.75"

    Returns a dict of float-valued parameters, or {} if parsing fails.
    """
    parts = config_label.split("_")

    if motion == "trnl":
        if len(parts) != 4:
            return {}
        _, x, y, z = parts
        return {"x": float(x), "y": float(y), "z": float(z)}

    if motion == "radial":
        # config label has the same n_x_y_z structure as trnl
        if len(parts) != 4:
            return {}
        _, x, y, z = parts
        return {"x": float(x), "y": float(y), "z": float(z)}

    if motion == "rot":
        if len(parts) != 5:
            return {}
        _, a, b, g, r = parts
        return {"alpha": float(a), "beta": float(b), "gamma": float(g), "R": float(r)}

    return {}


def missing_energies(E_scf, E_mp2, E_ccsdt):
    """
    Identify which energy components are missing from a parsed output file.

    Parameters
    ----------
    E_scf : float or None
        Total SCF energy extracted from the output file, or None if not found.
    E_mp2 : float or None
        Total MP2 energy extracted from the output file, or None if not found.
    E_ccsdt : float or None
        Total CCSD(T) energy extracted from the output file, or None if not found.

    Returns
    -------
    list of str
        Labels of missing energies, e.g. ["MP2", "CCSD(T)"] if only SCF
        converged. An empty list means all three energies were found.

    Notes
    -----
    A non-empty return value flags a non-converged or incomplete calculation.
    The corresponding configuration is excluded from raw_data entirely and
    its failure record is appended to raw_data_noconv.
    """
    missing = []
    if E_scf   is None: missing.append("SCF")
    if E_mp2   is None: missing.append("MP2")
    if E_ccsdt is None: missing.append("CCSD(T)")
    return missing


# ---------------------------------------------------------------------------
# Core builder
# ---------------------------------------------------------------------------

def build_raw_data(molecule, verbose=True):
    """
    Parse all NWChem .out files for *molecule* from its tarballs.

    Walks every *.tar.gz file found under the molecule directory using
    iter_out_texts(), extracts SCF, MP2, and CCSD(T) energies from each
    .out file, and assembles three dictionaries described below.

    Parameters
    ----------
    molecule : str
        Name of the molecule directory to process. Must match a subdirectory
        of bsse_db/data/ that contains *.tar.gz files, e.g. "Ne", "HF", "H2O".
    verbose : bool, optional
        If True (default), print a convergence summary to stdout after
        processing. Reports which motion types had failed configurations,
        or confirms that all configurations converged.

    Returns
    -------
    raw_data : dict
        Nested dictionary of converged energies with structure:
            raw_data[molecule][cluster][motion][config_id][method_key][basis][file_key]
        where:
            molecule   : str, e.g. "Ne"
            cluster    : str, cluster size e.g. "2" for dimer
            motion     : str, one of "trnl", "radial", "rot"
            config_id  : str, integer index e.g. "0", "1", ..., "267"
            method_key : str, one of "SCF", "MP2", "CCSD_T"
            basis      : str, e.g. "aug-cc-pvdz" or "aug-cc-pvtz"
            file_key   : str, one of "A_A", "AB_AB", "A_AB", "B_AB"
        Innermost values are floats in Hartree:
            SCF    -> total SCF energy
            MP2    -> MP2 correlation energy   (E_mp2 - E_scf)
            CCSD_T -> CCSD(T) component energy (E_ccsdt - E_mp2)

    raw_data_noconv : dict
        Records configurations excluded due to missing or non-converged
        energies. Structure:
            raw_data_noconv[molecule][cluster][motion][config_id]["failures"]
        where "failures" is a list of dicts, each containing:
            "missing"       : list of str, e.g. ["MP2", "CCSD(T)"]
            "method_dir"    : str
            "basis"         : str
            "file_key"      : str
            "tar_path"      : str
            "internal_path" : str

    parameters : dict
        Geometric parameters for each converged configuration. Structure:
            parameters[molecule][cluster][motion][config_id]
                -> {"config_label": str, "params": dict}
        config_label is the full folder name, e.g. "0_0.25_0.25_1.5".
        params is a dict of floats:
            trnl/radial : {"x": float, "y": float, "z": float}
            rot         : {"alpha": float, "beta": float,
                           "gamma": float, "R": float}
        Only populated for configurations present in raw_data.

    Exclusion rule
    --------------
    If ANY .out file in a (cluster, motion, config_id, method_dir, basis)
    folder is missing one or more of SCF/MP2/CCSD(T) energies, the entire
    folder is excluded from raw_data and logged in raw_data_noconv.
    """
    raw_data        = {}
    raw_data_noconv = {}
    parameters      = {}

    bucket              = {}   # basis_folder_id -> list of (file_key, E_scf, mp2_corr, ccsdt_comp)
    bad                 = set()
    config_label_by_cfg = {}   # (cluster, motion, config_key) -> config_label

    for tar_path, internal_path, text in iter_out_texts(molecule):
        internal_path = PurePosixPath(str(internal_path))

        if len(internal_path.parts) < 2:
            continue

        cluster      = internal_path.parts[0]
        motion       = motion_label(internal_path.parts[1])
        config_label = internal_path.parent.parent.parent.name
        config_key   = config_label.split("_", 1)[0]

        config_label_by_cfg[(cluster, motion, config_key)] = config_label

        file_key   = internal_path.stem
        basis      = internal_path.parent.name
        method_dir = internal_path.parent.parent.name

        basis_folder_id = (cluster, motion, config_key, method_dir, basis)

        E_scf   = extract_last(SCF_PATTERN,   text)
        E_mp2   = extract_last(MP2_PATTERN,   text)
        E_ccsdt = extract_last(CCSDT_PATTERN, text)

        if (E_scf is None) or (E_mp2 is None) or (E_ccsdt is None):
            bad.add(basis_folder_id)
            miss    = missing_energies(E_scf, E_mp2, E_ccsdt)
            cfg_bad = ensure_path(raw_data_noconv, (molecule, cluster, motion, config_key))
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

    # --- Commit only fully-converged folders ---
    converged_cfgs = set()

    for (cluster, motion, config_key, method_dir, basis), records in bucket.items():
        if (cluster, motion, config_key, method_dir, basis) in bad:
            continue

        converged_cfgs.add((cluster, motion, config_key))

        for file_key, E_scf, mp2_corr, ccsdt_component in records:
            cfg_node = ensure_path(raw_data, (molecule, cluster, motion, config_key))
            ensure_path(cfg_node, ("SCF",    basis))[file_key] = E_scf
            ensure_path(cfg_node, ("MP2",    basis))[file_key] = mp2_corr
            ensure_path(cfg_node, ("CCSD_T", basis))[file_key] = ccsdt_component

    # --- Build parameters for converged configs only ---
    for (cluster, motion, config_key) in converged_cfgs:
        config_label = config_label_by_cfg.get((cluster, motion, config_key))
        if config_label is None:
            continue
        node = ensure_path(parameters, (molecule, cluster, motion))
        node[config_key] = {
            "config_label": config_label,
            "params":       parse_params(motion, config_label),
        }

    # --- Verbose report ---
    if verbose:
        if raw_data_noconv:
            print(f"\n⚠️  WARNING: Some {molecule} configurations did NOT converge.")
            print("    Inspect `raw_data_noconv` for details.\n")
            for mol, mol_data in raw_data_noconv.items():
                for cluster, cluster_data in mol_data.items():
                    for motion, motion_data in cluster_data.items():
                        bad_cfgs = sorted(motion_data.keys(), key=int)
                        print(f"    {mol} | cluster {cluster} | {motion} | failed configs: {bad_cfgs}")
        else:
            print(f"\n✅  {molecule}: all configurations converged successfully.\n")

    return raw_data, raw_data_noconv, parameters


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Parse NWChem .out files and build raw_data pickle(s)."
    )
    parser.add_argument(
        "--mol",
        nargs="+",
        choices=["H2O", "Ne", "HF"],
        required=True,
        metavar="MOLECULE",
        help="Molecule(s) to process. Choices: H2O, Ne, HF. "
             "Examples: --mol Ne   --mol Ne HF   --mol H2O Ne HF",
    )
    args = parser.parse_args()

    for mol in args.mol:
        print(f"\n{'='*60}")
        print(f"  Processing: {mol}")
        print(f"{'='*60}")

        raw_data, raw_data_noconv, parameters = build_raw_data(mol, verbose=True)

        out_pickle = f"raw_data_{mol}.pickle"
        with open(out_pickle, "wb") as fh:
            pickle.dump(
                {"raw_data": raw_data, "raw_data_noconv": raw_data_noconv, "parameters": parameters},
                fh
            )
        print(f"  💾  Saved: {out_pickle}\n")