"""
two_body_bsse.py

Purpose
-------
Compute two-body BSSE for one or more molecules from their
many_body_interactions_{mol}.pickle files.

For each configuration, method, and basis set the BSSE is:

    BSSE = ΔE_IJ(I,J,IJ)  -  ΔE_IJ(IJ)
         = no_vmfc energy  -  vmfc energy

The result is stored in a nested dictionary that mirrors the
structure of the many_body_interactions pickle:

    two_body_bsse[mol][cluster][motion][config][method][basis]
        -> {'no_vmfc': float, 'vmfc': float, 'bsse': float}

Usage (CLI)
-----------
    python3 two_body_bsse.py --mol Ne
    python3 two_body_bsse.py --mol HF
    python3 two_body_bsse.py --mol Ne HF
    python3 two_body_bsse.py --mol H2O Ne HF

Output
------
    two_body_bsse_{mol}.pickle   (one file per molecule)

Import
------
    from two_body_bsse import twobody_bsse
    result = twobody_bsse(pickle_path, "Ne")
"""

import argparse
import pickle
from pathlib import Path

from open_pickle_data import open_pickle


# ---------------------------------------------------------------------------
# Core function (importable by dataframe-building scripts)
# ---------------------------------------------------------------------------

def twobody_bsse(pickle_file, molecule):
    """
    Compute two-body BSSE for a single molecule from its
    many_body_interactions pickle.

    Parameters
    ----------
    pickle_file : str or Path
        Path to many_body_interactions_{mol}.pickle
    molecule : str
        Molecule key, e.g. 'Ne', 'HF', 'H2O'

    Returns
    -------
    dict
        Nested dict:
        two_body_bsse[molecule][cluster][motion][config][method][basis]
            -> {'no_vmfc': float, 'vmfc': float, 'bsse': float}
    """
    raw_data = open_pickle(pickle_file)
    two_body_no_vmfc = raw_data["2b_no_vmfc"]
    two_body_vmfc    = raw_data["2b_vmfc"]

    two_body_bsse = {}
    two_body_bsse[molecule] = {}

    for cluster in two_body_no_vmfc[molecule]:
        two_body_bsse[molecule][cluster] = {}

        for motion in two_body_no_vmfc[molecule][cluster]:
            two_body_bsse[molecule][cluster][motion] = {}

            for config in two_body_no_vmfc[molecule][cluster][motion]:
                two_body_bsse[molecule][cluster][motion][config] = {}

                for method in two_body_no_vmfc[molecule][cluster][motion][config]:
                    two_body_bsse[molecule][cluster][motion][config][method] = {}

                    for basis in two_body_no_vmfc[molecule][cluster][motion][config][method]:
                        DE_no_vmfc = two_body_no_vmfc[molecule][cluster][motion][config][method][basis]["energy"]
                        DE_vmfc    = two_body_vmfc   [molecule][cluster][motion][config][method][basis]["energy"]

                        two_body_bsse[molecule][cluster][motion][config][method][basis] = {
                            "no_vmfc": DE_no_vmfc,
                            "vmfc":    DE_vmfc,
                            "bsse":    DE_no_vmfc - DE_vmfc,
                        }

    return two_body_bsse


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Compute two-body BSSE for one or more molecules."
    )
    parser.add_argument(
        "--mol",
        nargs="+",
        required=True,
        metavar="MOL",
        help="Molecule(s) to process, e.g. --mol Ne  or  --mol Ne HF H2O",
    )
    args = parser.parse_args()

    data_dir = Path(__file__).resolve().parent

    for mol in args.mol:
        pickle_in  = data_dir / f"many_body_interactions_{mol}.pickle"
        pickle_out = data_dir / f"two_body_bsse_{mol}.pickle"

        if not pickle_in.exists():
            print(f"[SKIP] {pickle_in.name} not found — skipping {mol}.")
            continue

        print(f"\n[{mol}] Reading {pickle_in.name} ...")
        result = twobody_bsse(pickle_in, mol)

        # Quick sanity print
        mol_data = result[mol]
        for cluster in mol_data:
            for motion in mol_data[cluster]:
                n_cfg = len(mol_data[cluster][motion])
                print(f"  cluster={cluster}  motion={motion:8s}  configs={n_cfg}")

        with open(pickle_out, "wb") as fh:
            pickle.dump(result, fh)
        print(f"💾  Saved: {pickle_out.name}")


if __name__ == "__main__":
    main()