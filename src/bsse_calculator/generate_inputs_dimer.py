from bsse_calculator.prepare_input import prepare_input


def dimer_geometry_str(geometry):
    ret = ""
    for molecule in geometry:
        for line in molecule:
            ret += f"  {line}\n"
    return ret


def dimer_input(geometry, input_file_path):
    """
    Creates NWChem input file for geometry to calculate E(AB, AB)
    """
    prepare_input(dimer_geometry_str(geometry), input_file_path)
