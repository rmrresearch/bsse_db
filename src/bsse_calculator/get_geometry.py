from math import sqrt


def get_geometry(HF_bond_length, **kwargs):
    if HF_bond_length < 0.5:
        raise ValueError("HF_bond_length must be greater than 0.5")
    if "FF_distance" in kwargs:
        FF_distance = kwargs["FF_distance"]
        if FF_distance < 1.5:
            raise ValueError("FF_distance must be greater than 1.5")
        return [
            [
                f"H 0 0 -{HF_bond_length}",
                "F 0 0 0",
            ],
            [f"F 0 0 {FF_distance}", f"H 0 0 {FF_distance + HF_bond_length}"],
        ]
    elif "F2_coords" in kwargs:
        F2_coords = kwargs["F2_coords"]
        if sqrt(F2_coords[0] ** 2 + F2_coords[1] ** 2 + F2_coords[2] ** 2) < 1.5:
            raise ValueError("F2_coords must be greater than 1.5")
        return [
            [f"H 0 0 -{HF_bond_length}", "F 0 0 0"],
            [
                f"F {F2_coords[0]} {F2_coords[1]} {F2_coords[2]}",
                f"H {F2_coords[0]} {F2_coords[1]} {(F2_coords[2] + HF_bond_length)}",
            ],
        ]
    else:
        raise ValueError("Please specify either FF_distance or F2_coords")
