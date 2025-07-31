from bsse_calculator.prepare_input import prepare_input
import numpy as np


def dimer_geometry_str(geometry):
    ret = ""
    for molecule in geometry:
        for line in molecule:
            ret += f"  {line}\n"
    return ret


def dimer_input(geometry, input_file_path, basis_set="aug-cc-pvdz"):
    """
    Creates NWChem input file for geometry to calculate E(AB, AB)
    """
    prepare_input(dimer_geometry_str(geometry), input_file_path, "energy", basis_set)
