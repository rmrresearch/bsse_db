import subprocess, os
from rdkit import Chem
import py3Dmol


def get_geometry(HF_bond_length, **kwargs):
    if "FF_distance" in kwargs:
        FF_distance = kwargs["FF_distance"]
        return [
            [
                f"H 0 0 -{HF_bond_length}",
                "F 0 0 0",
            ],
            [f"F 0 0 {FF_distance}", f"H 0 0 {FF_distance + HF_bond_length}"],
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


def dimer_input(geometry, input_file_path):
    """
    Creates nwchem input file for geometry to calculate E(AB, AB)
    """
    prepare_input(dimer_geometry_str(geometry), input_file_path)


def bsse_corrected_monomer_input(geometry, input_A, input_B):
    """
    Creates nwchem input files for geometry to calculate E(A, AB) and E(B, AB)
    """
    prepare_input(bsse_corrected_monomer_geometry_str(geometry, 0), input_A)
    prepare_input(bsse_corrected_monomer_geometry_str(geometry, 1), input_B)


def monomer_input(geometry, input_file_path):
    """Creates nwchem input file for geometry to calculate E(A, A)"""
    prepare_input(monomer_geometry_str(geometry), input_file_path)


def get_bsse(geometry, scripts_dir="scripts", force_rerun=True) -> dict:
    """
    Calculates BSSE for geometry (of dimer)
    Returns a dictionary of the form
    {
    "Delta_E_AB_AB": float,
    "BSSE_A": float,
    "BSSE_B": float
    }
    where Delta_E_AB_AB = ΔE(AB, AB) and BSSE_A = ε(A, AB) and BSSE_B = ε(B, AB)
    """
    os.makedirs(scripts_dir, exist_ok=True)
    dimer_input(geometry, os.path.join(scripts_dir, "input_dimer.txt"))
    bsse_corrected_monomer_input(
        geometry,
        os.path.join(scripts_dir, "input_A_AB.txt"),
        os.path.join(scripts_dir, "input_B_AB.txt"),
    )
    monomer_input(geometry, os.path.join(scripts_dir, "input_monomer.txt"))
    if force_rerun:
        run_software(
            [
                os.path.join(scripts_dir, "input_dimer.txt"),
                os.path.join(scripts_dir, "input_A_AB.txt"),
                os.path.join(scripts_dir, "input_B_AB.txt"),
                os.path.join(scripts_dir, "input_monomer.txt"),
            ]
        )
    else:
        needed_files = [
            os.path.join(scripts_dir, "input_dimer.txt"),
            os.path.join(scripts_dir, "input_A_AB.txt"),
            os.path.join(scripts_dir, "input_B_AB.txt"),
            os.path.join(scripts_dir, "input_monomer.txt"),
        ]
        if os.path.isdir(scripts_dir) and set(needed_files).issubset(
            [os.path.join(scripts_dir, file) for file in os.listdir(scripts_dir)]
        ):
            pass
        else:
            run_software(
                [
                    os.path.join(scripts_dir, "input_dimer.txt"),
                    os.path.join(scripts_dir, "input_A_AB.txt"),
                    os.path.join(scripts_dir, "input_B_AB.txt"),
                    os.path.join(scripts_dir, "input_monomer.txt"),
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
    return {
        "Delta_E_AB_AB": total_energy - E_A_AB - E_B_AB,
        "BSSE_A": E_monomer - E_A_AB,
        "BSSE_B": E_monomer - E_B_AB,
    }


def prepare_input(geometry_str, input_file_path):
    with open(input_file_path, "w") as f:
        f.write("geometry\n")
        f.write(geometry_str)
        f.write("end\n")
        f.write("basis spherical\n")
        f.write("  F library aug-cc-pvdz\n")
        f.write("  H library aug-cc-pvdz\n")
        f.write("end\n")
        f.write("task ccsd(t) energy\n")


def dimer_geometry_str(geometry):
    ret = ""
    for molecule in geometry:
        for line in molecule:
            ret += f"  {line}\n"
    return ret


def bsse_corrected_monomer_geometry_str(geometry, ghost_monomer):
    ret = ""
    for i, molecule in enumerate(geometry):
        for line in molecule:
            ret += f"  {line}\n" if i == ghost_monomer else f"  bq{line}\n"
    return ret


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
            if "Total CCSD(T) energy" in line:
                total_energy = float(line.split()[-1])
            if "(T) corr. energy" in line:
                paren_t_correlation_energy = float(line.split()[-1])
            if "CCSD corr. energy" in line:
                CSSD_corr_energy = float(line.split()[-1])
            if "MP2 Corr. energy:" in line:
                mp2_correlation_energy = float(line.split()[-1])
            if "Total SCF energy" in line:
                scf_energy = float(line.split()[-1])
        return {
            "Total_energy": total_energy,
            "MP2_correlation_energy": mp2_correlation_energy,
            "(T)_correlation_energy": paren_t_correlation_energy,
            "CCSD corr.energy": CSSD_corr_energy,
            "SCF_energy": scf_energy,
        }


def monomer_geometry_str(geometry):
    monomer = geometry[0]
    ret = ""
    for line in monomer:
        ret += f"  {line}\n"
    return ret


def run_software(input_files, force_run=True):
    """
    Runs nwchem on the input files.

    Args:
        input_files (list): list of input files
        force_run (bool): whether to run nwchem or not if output files already exist
    Returns the directory of the output path"""
    if force_run:
        for input_file in input_files:
            output_file = input_file.replace("input", "output")
            with open(output_file, "w") as out:
                subprocess.run(["nwchem", input_file], stdout=out)
    else:
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
                output_file = input_file.replace("input", "output")
                with open(output_file, "w") as out:
                    subprocess.run(["nwchem", input_file], stdout=out)


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
