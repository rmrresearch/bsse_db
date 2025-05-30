import subprocess
import os


def get_geometry(HF_bond_length, **kwargs):
    if "FF_distance" in kwargs:
        FF_distance = kwargs["FF_distance"]
        return [
            f"H 0.0 0.0 -{HF_bond_length}",
            "F 0.0 0.0 0.0",
            f"F 0.0 0.0 {FF_distance}",
            f"H 0.0 0.0 {FF_distance + HF_bond_length}",
        ]
    elif "F2_coords" in kwargs:
        F2_coords = kwargs["F2_coords"]
        return [
            [f"H 0.0 0.0 -{HF_bond_length}", "F 0.0 0.0 0.0"],
            [
                f"F {F2_coords[0]:.1f} {F2_coords[1]:.1f} {F2_coords[2]:.1f}",
                f"H {F2_coords[0]:.1f} {F2_coords[1]:.1f} {(F2_coords[2] + HF_bond_length):.1f}",
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
    prepare_input(dimer_geometry_str(geometry), input_file_path)


def bsse_corrected_monomer_input(geometry, input_A, input_B):
    prepare_input(bsse_corrected_monomer_geometry_str(geometry, 0), input_A)
    prepare_input(bsse_corrected_monomer_geometry_str(geometry, 1), input_B)


def read_total_energy(file_path):
    with open(file_path, "r") as f:
        for line in f:
            if "Total Energy =" in line:
                return float(line.split()[-1])
        raise ValueError(f"Total energy not found in {file_path}")


def monomer_geometry_str(geometry):
    monomer = geometry[0]
    ret = ""
    for line in monomer:
        ret += f"  {line}\n"
    return ret


def monomer_input(geometry, input_file_path):
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
    Returns a dictionary of the form
    {
      "Delta_E_AB_AB": float,
      "BSSE_A": float,
      "BSSE_B": float
    }"""
    prev_dir = os.getcwd()
    os.makedirs(scripts_dir, exist_ok=True)
    os.chdir(scripts_dir)
    dimer_input(geometry, "input_dimer.txt")
    bsse_corrected_monomer_input(geometry, "input_A_AB.txt", "input_B_AB.txt")
    monomer_input(geometry, "input_monomer.txt")
    output_path = run_psi4(
        ["input_dimer.txt", "input_A_AB.txt", "input_B_AB.txt", "input_monomer.txt"]
    )
    total_energy = read_total_energy(os.path.join(output_path, "output_dimer.txt"))
    E_A_AB = read_total_energy(os.path.join(output_path, "output_A_AB.txt"))
    E_B_AB = read_total_energy(os.path.join(output_path, "output_B_AB.txt"))
    E_monomer = read_total_energy(os.path.join(output_path, "output_monomer.txt"))
    os.chdir(prev_dir)
    return {
        "Delta_E_AB_AB": total_energy - E_A_AB - E_B_AB,
        "BSSE_A": E_monomer - E_A_AB,
        "BSSE_B": E_monomer - E_B_AB,
    }
