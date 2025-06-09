from bsse_calculator.prepare_input import prepare_input


def monomer_geometry_str(geometry):
    monomer = geometry[0]
    ret = ""
    for line in monomer:
        ret += f"  {line}\n"
    return ret


def monomer_input(geometry, input_file_path):
    """Creates NWChem input file for geometry to calculate E(A, A)"""
    prepare_input(monomer_geometry_str(geometry), input_file_path)
