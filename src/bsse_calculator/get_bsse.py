import os, subprocess
from bsse_calculator.generate_inputs_monomer import monomer_input
from bsse_calculator.generate_inputs_dimer import dimer_input
from bsse_calculator.generate_bsse_corrected_monomer_input import (
    bsse_corrected_monomer_input,
)
from bsse_calculator.read_total_energy import read_energies
import numpy as np
import tempfile
import shutil
import subprocess
from pathlib import Path
import glob


def run_software(input_files, force_run=True):
    """
    Runs NWChem on the input files.

    Args:
        input_files (list): list of input files
        force_run (bool): whether to run NWChem or not if output files already exist
    Returns the directory of the output path"""
    for input_path in input_files:
        input_path = Path(input_path).resolve()
        output_path = input_path.with_name(input_path.name.replace("input", "output"))
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir_path = Path(tmpdir)

            tmp_input = tmpdir_path / input_path.name
            shutil.copy(input_path, tmp_input)
            with open(output_path, "w") as out:
                result = subprocess.run(
                    ["mpirun", "--use-hwthread-cpus", "nwchem", tmp_input],
                    stdout=out,  # streamed to file
                    stderr=subprocess.PIPE,  # captured in memory
                    text=True,
                )
            if result.stderr:
                return
        for file_path in glob.glob("input_*"):
            os.remove(file_path)


def get_bsses(geometry, basis_set, scripts_dir="scripts", force_rerun=True) -> dict:
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
    dimer_input(geometry, os.path.join(scripts_dir, "input_dimer.txt"), basis_set)
    bsse_corrected_monomer_input(
        geometry,
        os.path.join(scripts_dir, "input_A_AB.txt"),
        os.path.join(scripts_dir, "input_B_AB.txt"),
        basis_set,
    )
    monomer_input(geometry, os.path.join(scripts_dir, "input_monomer.txt"), basis_set)
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
    try:
        dimer_energies = read_energies(os.path.join(scripts_dir, "output_dimer.txt"))
        E_A_AB_energies = read_energies(os.path.join(scripts_dir, "output_A_AB.txt"))
        E_B_AB_energies = read_energies(os.path.join(scripts_dir, "output_B_AB.txt"))
        E_monomer_energies = read_energies(
            os.path.join(scripts_dir, "output_monomer.txt")
        )
        output = {}
        for energy in dimer_energies.keys():
            output[f"{energy}_Delta_E_AB_AB"] = (
                dimer_energies[energy]
                - E_A_AB_energies[energy]
                - E_B_AB_energies[energy]
            )

            output[f"{energy}_BSSE_A"] = (
                E_monomer_energies[energy] - E_A_AB_energies[energy]
            )

            output[f"{energy}_BSSE_B"] = (
                E_monomer_energies[energy] - E_B_AB_energies[energy]
            )

        return output
    except:
        raise ValueError("Failed to calculate BSSE")
