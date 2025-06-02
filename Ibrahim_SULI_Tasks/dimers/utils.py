import subprocess
import os
from rdkit import Chem
import py3Dmol


def get_geometry(HF_bond_length, **kwargs):
    if "FF_distance" in kwargs:
        FF_distance = kwargs["FF_distance"]
        return [
            f"H 0 0 -{HF_bond_length}",
            "F 0 0 0",
            f"F 0 0 {FF_distance}",
            f"H 0 0 {FF_distance + HF_bond_length}",
        ]
    elif "F2_coords" in kwargs:
        F2_coords = kwargs["F2_coords"]
        return [
            [f"H 0 0 -{HF_bond_length}", "F 0 0 0"],
            [
                f"F {F2_coords[0]} {F2_coords[1]} {F2_coords[2]}",
                f"H {F2_coords[0]} {F2_coords[1]} {(F2_coords[2] + HF_bond_length)}",
            ],
        ]
    else:
        raise ValueError("Please specify either FF_distance or F2_coords")


def prepare_input(geometry_str, input_file_path):
    with open(input_file_path, "w") as f:
        f.write("memory 600 mb\n")
        f.write("\n")
        f.write("molecule HF {\n")
        f.write(geometry_str)
        f.write("}\n")
        f.write("\n")
        f.write("set basis aug-cc-pVDZ\n")
        f.write("energy('ccsd(t)')\n")


def dimer_geometry_str(geometry):
    ret = ""
    for i, molecule in enumerate(geometry):
        for line in molecule:
            ret += f"  {line}\n"
        if i != len(geometry) - 1:
            ret += "  --\n"
    return ret


def bsse_corrected_monomer_geometry_str(geometry, ghost_monomer):
    ret = ""
    for i, molecule in enumerate(geometry):
        for line in molecule:
            ret += f"  {line}\n" if i == ghost_monomer else f"  @{line}\n"
        if i != len(geometry) - 1:
            ret += f"  --\n"
    return ret


def dimer_input(geometry, input_file_path):
    """
    Creates psi4 input file for geometry to calculate E(AB, AB)
    """
    prepare_input(dimer_geometry_str(geometry), input_file_path)


def bsse_corrected_monomer_input(geometry, input_A, input_B):
    """
    Creates psi4 input files for geometry to calculate E(A, AB) and E(B, AB)
    """
    prepare_input(bsse_corrected_monomer_geometry_str(geometry, 0), input_A)
    prepare_input(bsse_corrected_monomer_geometry_str(geometry, 1), input_B)


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


def monomer_geometry_str(geometry):
    monomer = geometry[0]
    ret = ""
    for line in monomer:
        ret += f"  {line}\n"
    return ret


def monomer_input(geometry, input_file_path):
    """Creates psi4 input file for geometry to calculate E(A, A)"""
    prepare_input(monomer_geometry_str(geometry), input_file_path)


def run_psi4(input_files):
    """
    Runs psi4 on the input files.
    Returns the directory of the output path"""
    for input_file in input_files:
        subprocess.run(["psi4", input_file, input_file.replace("input", "output")])
    return (
        os.getcwd()
    )  # this is done so we can mock the method to return the test asserts folder for testing


def get_bsse(geometry, scripts_dir="scripts") -> dict:
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
    output_path = run_psi4(
        ["input_dimer.txt", "input_A_AB.txt", "input_B_AB.txt", "input_monomer.txt"]
    )
    total_energy = read_energies(os.path.join(output_path, "output_dimer.txt"))[
        "Total_energy"
    ]
    E_A_AB = read_energies(os.path.join(output_path, "output_A_AB.txt"))["Total_energy"]
    E_B_AB = read_energies(os.path.join(output_path, "output_B_AB.txt"))["Total_energy"]
    E_monomer = read_energies(os.path.join(output_path, "output_monomer.txt"))[
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
