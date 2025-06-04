import subprocess
import os
from rdkit import Chem
import py3Dmol
from bsse_calculator.generate_inputs import (
    dimer_input,
    monomer_input,
    bsse_corrected_monomer_input,
)


def read_energies(file_path):
    """
    Returns a dictionary of the form
    {
      "Total_energy": float,
      "MP2_correlation_energy": float,
      "CCSD(T)_correlation_energy": float
    }
    output
    """
    with open(file_path, "r") as f:
        for line in f:
            if "Total Energy =" in line:
                total_energy = float(line.split()[-1])
            if "MP2 correlation energy" in line:
                mp2_correlation_energy = float(line.split()[-1])
            if "CCSD(T) total energy" in line:
                ccsd_total_energy = float(line.split()[-1])
            if "SCF energy" in line:
                scf_energy = float(line.split()[-1])
        return {
            "Total_energy": total_energy,
            "MP2_correlation_energy": mp2_correlation_energy,
            "CCSD(T)_correlation_energy": ccsd_total_energy
            - scf_energy
            - mp2_correlation_energy,
            "SCF_energy": scf_energy,
        }


def run_psi4(input_files, force_run=True):
    """
    Runs psi4 on the input files.

    Args:
        input_files (list): list of input files
        force_run (bool): whether to run psi4 or not if output files already exist
    Returns the directory of the output path"""
    if not force_run:
        needed_files = [
            "output_dimer.txt",
            "output_A_AB.txt",
            "output_B_AB.txt",
            "output_monomer.txt",
        ]
        if set(needed_files).issubset(os.listdir(os.getcwd())):
            return
        else:
            for input_file in input_files:
                subprocess.run(
                    ["psi4", input_file, input_file.replace("input", "output")]
                )


def get_bsse(geometry, scripts_dir="scripts", force_rerun=True) -> dict:
    """
    Calculates BSSE for geometry (of dimer)
    Returns a dictonary of the form
    {
      "Delta_E_AB_AB": float,
      "BSSE_A": float,
      "BSSE_B": float
    }
    where Delta_E_AB_AB = ΔE(AB, AB) and BSSE_A = ε(A, AB) and BSSE_B = ε(B, AB)
    """
    prev_dir = os.getcwd()
    os.makedirs(scripts_dir, exist_ok=True)
    os.chdir(scripts_dir)
    dimer_input(geometry, "input_dimer.txt")
    bsse_corrected_monomer_input(geometry, "input_A_AB.txt", "input_B_AB.txt")
    monomer_input(geometry, "input_monomer.txt")
    if force_rerun:
        run_psi4(
            ["input_dimer.txt", "input_A_AB.txt", "input_B_AB.txt", "input_monomer.txt"]
        )
    else:
        needed_files = [
            os.path.join(scripts_dir, file)
            for file in [
                "output_dimer.txt",
                "output_A_AB.txt",
                "output_B_AB.txt",
                "output_monomer.txt",
            ]
        ]
        if os.path.isdir(scripts_dir) and set(needed_files).issubset(
            os.listdir(scripts_dir)
        ):
            output_path = scripts_dir
        else:
            run_psi4(
                [
                    "input_dimer.txt",
                    "input_A_AB.txt",
                    "input_B_AB.txt",
                    "input_monomer.txt",
                ]
            )
    total_energy = read_energies(os.path.join(scripts_dir, "output_dimer.txt"))[
        "Total_energy"
    ]
    E_A_AB = read_energies(os.path.join(scripts_dir, "output_A_AB.txt"))["Total_energy"]
    E_B_AB = read_energies(os.path.join(scripts_dir, "output_B_AB.txt"))["Total_energy"]
    E_monomer = read_energies(os.path.join(scripts_dir, "output_monomer.txt"))[
        "Total_energy"
    ]
    os.chdir(prev_dir)
    return {
        "Delta_E_AB_AB": total_energy - E_A_AB - E_B_AB,
        "BSSE_A": E_monomer - E_A_AB,
        "BSSE_B": E_monomer - E_B_AB,
    }


def make_diagram(geometry, output_path):
    """
    Creates 3D diagram of geometry saved as html file
    """
    bonds = [
        (0, 1, Chem.BondType.SINGLE),
        (2, 3, Chem.BondType.SINGLE),
    ]
    mol = Chem.RWMol()
    for row in [*geometry[0], *geometry[1]]:
        atom_symbol = row.split()[0]
        mol.AddAtom(Chem.Atom(atom_symbol))
    for i, j, bond_type in bonds:
        mol.AddBond(i, j, bond_type)

    # Add conformer with coordinates
    conf = Chem.Conformer(mol.GetNumAtoms())
    for i, row in enumerate([*geometry[0], *geometry[1]]):
        _, x, y, z = row.split()
        conf.SetAtomPosition(i, Chem.rdGeometry.Point3D(float(x), float(y), float(z)))
    mol.AddConformer(conf)

    mol_block = Chem.MolToMolBlock(mol)

    view = py3Dmol.view(width=400, height=400)
    view.addModel(mol_block, "mol")
    view.setStyle({"stick": {}, "sphere": {"scale": 0.3}})
    view.zoomTo()
    html_str = view._make_html()
    with open(output_path, "w") as f:
        f.write(html_str)
