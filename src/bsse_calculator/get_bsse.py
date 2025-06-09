import os, subprocess
from bsse_calculator.generate_inputs_monomer import monomer_input
from bsse_calculator.generate_inputs_dimer import dimer_input
from bsse_calculator.generate_bsse_corrected_monomer_input import (
    bsse_corrected_monomer_input,
)
from bsse_calculator.read_total_energy import read_energies


def run_software(input_files, force_run=True):
    """
    Runs NWChem on the input files.

    Args:
        input_files (list): list of input files
        force_run (bool): whether to run NWChem or not if output files already exist
    Returns the directory of the output path"""
    for input_file in input_files:
        output_file = input_file.replace("input", "output")
        with open(output_file, "w") as out:
            subprocess.run(["nwchem", input_file], stdout=out)


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
            os.path.join(scripts_dir, "output_dimer.txt"),
            os.path.join(scripts_dir, "output_A_AB.txt"),
            os.path.join(scripts_dir, "input_B_AB.txt"),
            os.path.join(scripts_dir, "output_monomer.txt"),
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
