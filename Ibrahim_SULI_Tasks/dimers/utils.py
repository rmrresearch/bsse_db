import subprocess
import os


def get_geometry(HF_bond_length, FF_distance):
    return [
        [f"H 0.0 0.0 -{HF_bond_length}", "F 0.0 0.0 0.0"],
        [f"F 0.0 0.0 {FF_distance}", f"H 0.0 0.0 {FF_distance + HF_bond_length}"],
    ]


def dimer_input(geometry, input_file_path):
    with open(input_file_path, "w") as f:
        f.write("memory 600 mb\n")
        f.write("\n")
        f.write("molecule HF {\n")
        for i, molecule in enumerate(geometry):
            for line in molecule:
                f.write(f"  {line}\n")
            if i != len(geometry) - 1:
                f.write("  --\n")
        f.write("}\n")
        f.write("\n")
        f.write("set basis aug-cc-pVDZ\n")
        f.write("energy('ccsd(t)')\n")


def monomer_input(geometry, input_A, input_B):
    with open(input_A, "w") as f:
        f.write("memory 600 mb\n")
        f.write("\n")
        f.write("molecule HF {\n")
        for i, molecule in enumerate(geometry):
            for line in molecule:
                f.write(f"  {line}\n" if i == 0 else f"  @{line}\n")
            if i != len(geometry) - 1:
                f.write("  --\n")
        f.write("}\n")
        f.write("\n")
        f.write("set basis aug-cc-pVDZ\n")
        f.write("energy('ccsd(t)')\n")
    with open(input_B, "w") as f:
        f.write("memory 600 mb\n")
        f.write("\n")
        f.write("molecule HF {\n")
        for i, molecule in enumerate(geometry):
            for line in molecule:
                f.write(f"  {line}\n" if i == 1 else f"  @{line}\n")
            if i != len(geometry) - 1:
                f.write("  --\n")
        f.write("}\n")
        f.write("\n")
        f.write("set basis aug-cc-pVDZ\n")
        f.write("energy('ccsd(t)')\n")


def read_total_energy(file_path):
    with open(file_path, "r") as f:
        for line in f:
            if "Total Energy =" in line:
                return float(line.split()[-1])
        raise ValueError(f"Total energy not found in {file_path}")


def get_optimized_monomer_energy(geometry, input_file_path):
    monomer = geometry[0]
    with open(input_file_path, "w") as f:
        f.write("memory 600 mb\n")
        f.write("\n")
        f.write("molecule HF {\n")
        for line in monomer:
            f.write(f"  {line}\n")
        f.write("}\n")
        f.write("\n")
        f.write("set basis aug-cc-pVDZ\n")
        f.write("energy('ccsd(t)')\n")


def get_bsse(geometry, scripts_dir="scripts") -> dict:
    """
    Returns a dictonary of the form
    {
      "Delta_E_AB_AB": float,
      "BSSE_A": float,
      "BSSE_B": float
    }"""
    prev_dir = os.getcwd()
    os.makedirs(scripts_dir, exist_ok=True)
    os.chdir(scripts_dir)
    dimer_input(geometry, "input_dimer.txt")
    monomer_input(geometry, "input_A_AB.txt", "input_B_AB.txt")
    get_optimized_monomer_energy(geometry, "input_monomer.txt")
    subprocess.run(["psi4", "input_dimer.txt", "output_dimer.txt"])
    subprocess.run(["psi4", "input_A_AB.txt", "output_A_AB.txt"])
    subprocess.run(["psi4", "input_B_AB.txt", "output_B_AB.txt"])
    subprocess.run(["psi4", "input_monomer.txt", "output_monomer.txt"])
    total_energy = read_total_energy("output_dimer.txt")
    E_A_AB = read_total_energy("output_A_AB.txt")
    E_B_AB = read_total_energy("output_B_AB.txt")
    E_monomer = read_total_energy("output_monomer.txt")
    os.chdir(prev_dir)
    return {
        "Delta_E_AB_AB": total_energy - E_A_AB - E_B_AB,
        "BSSE_A": E_monomer - E_A_AB,
        "BSSE_B": E_monomer - E_B_AB,
    }
