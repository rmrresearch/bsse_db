from bsse_calculator.prepare_input import prepare_input


def bsse_corrected_monomer_geometry_str(geometry, ghost_monomer):
    ret = ""
    for i, molecule in enumerate(geometry):
        for line in molecule:
            ret += f"  {line}\n" if i == ghost_monomer else f"  bq{line}\n"
    return ret


def bsse_corrected_monomer_input(geometry, input_A, input_B, basis_set):
    """
    Creates NWChem input files for geometry to calculate E(A, AB) and E(B, AB)
    """
    prepare_input(
        bsse_corrected_monomer_geometry_str(geometry, 0), input_A, "energy", basis_set
    )
    prepare_input(
        bsse_corrected_monomer_geometry_str(geometry, 1), input_B, "energy", basis_set
    )
