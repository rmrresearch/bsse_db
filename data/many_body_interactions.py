"""
many_body_interactions.py
=========================
Constructs two-body interaction energy dictionaries for BSSE analysis
for any supported molecule (H2O, Ne, HF) via a command-line interface.

Purpose
-------
This script implements part of Task 4 of the PI task list: given the parsed
raw energies from build_raw_data.py, compute the two-body interaction energies
in both the standard (no-VMFC) and counterpoise-corrected (VMFC) formulations.
The results are serialized to a molecule-specific pickle file for downstream
use by the BSSE computation and dataframe-building scripts (Tasks 5 and 6).

Two interaction energy definitions
------------------------------------
1. Standard (no-VMFC) two-body interaction energy:

       ΔE_IJ(I,J,IJ) = E_IJ(IJ) - E_I(I) - E_J(J)

   Monomer energies are computed in their own (monomer) basis sets.

2. Counterpoise / VMFC two-body interaction energy:

       ΔE_IJ(IJ) = E_IJ(IJ) - E_I(IJ) - E_J(IJ)

   Monomer energies are computed in the full dimer basis set (with ghost
   functions on the other monomer).

Supported molecules
-------------------
H2O  — dimer: translational (trnl)
Ne   — dimer: translational (trnl) + radial
HF   — dimer: translational (trnl) + radial + rotational (rot)

Note: all motion types (trnl, radial, rot) are handled dynamically by the
nested loops — no motion-specific branching is required.

Homodimer / identical monomer handling
---------------------------------------
For Ne and HF dimers, monomers A and B have identical geometries, so no
separate B_B calculation was performed. When B_B is absent from the raw
data, E_B(B) is set equal to E_A(A). This is correct by symmetry.

Input
-----
raw_data_{MOL}.pickle — produced by build_raw_data.py
    Expected structure:
        data["raw_data"][molecule][cluster][motion][config][method][basis][file_key]

Output
------
many_body_interactions_{MOL}.pickle — one file per molecule
    Contains a dict with two keys:
        "2b_no_vmfc" : nested dict of standard interaction energies
        "2b_vmfc"    : nested dict of VMFC interaction energies

    Both dicts share the structure:
        result[molecule][cluster][motion][config][method][basis]["energy"] = float

Output dict structure
---------------------
result[molecule][cluster][motion][config][method][basis]["energy"]
    where:
        molecule : str, e.g. "Ne" or "HF"
        cluster  : str, e.g. "2" for dimer
        motion   : str, one of "trnl", "radial", "rot"
        config   : str, integer config index e.g. "0", "1", ...
        method   : str, one of "SCF", "MP2", "CCSD_T"
        basis    : str, e.g. "aug-cc-pvdz" or "aug-cc-pvtz"
        "energy" : float, interaction energy in Hartree

Usage
-----
Run from bsse_db/data/ with the raw_data pickle(s) present:

    # single molecule
    python3 many_body_interactions.py --mol Ne

    # two molecules
    python3 many_body_interactions.py --mol Ne HF

    # all three
    python3 many_body_interactions.py --mol H2O Ne HF

Dependencies
------------
    open_pickle_data.py — must be in the same directory; provides
                          open_pickle(path) which loads and returns
                          the contents of a pickle file.
"""

import pickle
import argparse
from pathlib import Path
from open_pickle_data import open_pickle


# ---------------------------------------------------------------------------
# Two-body interaction energy builders
# ---------------------------------------------------------------------------

def _dimer_no_vmfc(pickle_file, molecule):
    """
    Construct standard (no-VMFC) two-body interaction energies for a dimer.

    Computes:
        ΔE_IJ(I,J,IJ) = E_IJ(IJ) - E_I(I) - E_J(J)

    For homodimers or dimers with identical monomer geometries (Ne, HF),
    B_B is absent from the raw data and E_B(B) is set equal to E_A(A).

    Parameters
    ----------
    pickle_file : str or Path
        Path to raw_data_{mol}.pickle produced by build_raw_data.py.
    molecule : str
        Molecule label, e.g. "Ne", "HF", "H2O".

    Returns
    -------
    dict
        Nested dict of standard interaction energies:
            result[molecule][cluster][motion][config][method][basis]["energy"]
    """
    raw_data = open_pickle(pickle_file)['raw_data']
    two_body_no_vmfc = {molecule: {}}

    for cluster in raw_data[molecule]:
        two_body_no_vmfc[molecule][cluster] = {}

        for motion in raw_data[molecule][cluster]:
            two_body_no_vmfc[molecule][cluster][motion] = {}

            for config in raw_data[molecule][cluster][motion]:
                two_body_no_vmfc[molecule][cluster][motion][config] = {}

                for method in raw_data[molecule][cluster][motion][config]:
                    two_body_no_vmfc[molecule][cluster][motion][config][method] = {}

                    for basis in raw_data[molecule][cluster][motion][config][method]:
                        two_body_no_vmfc[molecule][cluster][motion][config][method][basis] = {}

                        file_data = raw_data[molecule][cluster][motion][config][method][basis]
                        E_A_A   = file_data['A_A']
                        E_AB_AB = file_data['AB_AB']

                        # B_B absent for Ne and HF (identical monomer geometries)
                        E_B_B = file_data.get('B_B', E_A_A)

                        energy = E_AB_AB - (E_A_A + E_B_B)
                        two_body_no_vmfc[molecule][cluster][motion][config][method][basis]['energy'] = energy

    return two_body_no_vmfc


def _dimer_vmfc(pickle_file, molecule):
    """
    Construct counterpoise / VMFC two-body interaction energies for a dimer.

    Computes:
        ΔE_IJ(IJ) = E_IJ(IJ) - E_I(IJ) - E_J(IJ)

    Monomer energies are computed in the full dimer basis (ghost functions
    on the other monomer), providing the Boys-Bernardi counterpoise correction.

    Parameters
    ----------
    pickle_file : str or Path
        Path to raw_data_{mol}.pickle produced by build_raw_data.py.
    molecule : str
        Molecule label, e.g. "Ne", "HF", "H2O".

    Returns
    -------
    dict
        Nested dict of VMFC interaction energies:
            result[molecule][cluster][motion][config][method][basis]["energy"]
    """
    raw_data = open_pickle(pickle_file)['raw_data']
    two_body_vmfc = {molecule: {}}

    for cluster in raw_data[molecule]:
        two_body_vmfc[molecule][cluster] = {}

        for motion in raw_data[molecule][cluster]:
            two_body_vmfc[molecule][cluster][motion] = {}

            for config in raw_data[molecule][cluster][motion]:
                two_body_vmfc[molecule][cluster][motion][config] = {}

                for method in raw_data[molecule][cluster][motion][config]:
                    two_body_vmfc[molecule][cluster][motion][config][method] = {}

                    for basis in raw_data[molecule][cluster][motion][config][method]:
                        two_body_vmfc[molecule][cluster][motion][config][method][basis] = {}

                        file_data = raw_data[molecule][cluster][motion][config][method][basis]
                        E_A_AB  = file_data['A_AB']
                        E_B_AB  = file_data['B_AB']
                        E_AB_AB = file_data['AB_AB']

                        energy = E_AB_AB - (E_A_AB + E_B_AB)
                        two_body_vmfc[molecule][cluster][motion][config][method][basis]['energy'] = energy

    return two_body_vmfc


# ---------------------------------------------------------------------------
# Main builder and serializer
# ---------------------------------------------------------------------------

def many_body_int_pickle(pickle_file, molecule, outdir="."):
    """
    Compute and serialize two-body interaction energies for one molecule.

    Builds both the standard (no-VMFC) and counterpoise (VMFC) two-body
    interaction energy dictionaries and writes them to a molecule-specific
    pickle file named:

        many_body_interactions_{molecule}.pickle

    Parameters
    ----------
    pickle_file : str or Path
        Path to raw_data_{mol}.pickle produced by build_raw_data.py.
    molecule : str
        Molecule label, e.g. "Ne", "HF", "H2O".
    outdir : str or Path, optional
        Directory where the output pickle will be written.
        Defaults to the current directory.

    Output
    ------
    Writes many_body_interactions_{molecule}.pickle containing:
        {
            "2b_no_vmfc": {molecule: ...},
            "2b_vmfc":    {molecule: ...}
        }

    Notes
    -----
    - Overwrites any existing file with the same name.
    - All motion types (trnl, radial, rot) are handled automatically
      by the dynamic nested loops — no motion-specific branching needed.
    - For Ne and HF, the absence of B_B in the raw data is handled
      inside _dimer_no_vmfc() by falling back to E_A_A.
    """
    two_body_no_vmfc = _dimer_no_vmfc(pickle_file, molecule)
    two_body_vmfc    = _dimer_vmfc(pickle_file, molecule)

    body_interactions = {
        "2b_no_vmfc": two_body_no_vmfc,
        "2b_vmfc":    two_body_vmfc,
    }

    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    out_file = outdir / f"many_body_interactions_{molecule}.pickle"
    with out_file.open("wb") as f:
        pickle.dump(body_interactions, f)

    print(f"  💾  Saved: {out_file}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Compute two-body interaction energies (no-VMFC and VMFC) "
                    "and serialize to molecule-specific pickle files."
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

        pickle_file = f"raw_data_{mol}.pickle"
        many_body_int_pickle(pickle_file, mol, outdir=".")