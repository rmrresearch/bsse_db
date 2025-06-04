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


def monomer_geometry_str(geometry):
    monomer = geometry[0]
    ret = ""
    for line in monomer:
        ret += f"  {line}\n"
    return ret


def monomer_input(geometry, input_file_path):
    """Creates psi4 input file for geometry to calculate E(A, A)"""
    prepare_input(monomer_geometry_str(geometry), input_file_path)
