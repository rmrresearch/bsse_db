from math import sqrt
import numpy as np
from bsse_calculator.generate_nmer_geometry import check_no_overlaps


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


def euler_angle_rotation(phi, theta, psi):
    rotation_z = np.array(
        [[np.cos(psi), -np.sin(psi), 0], [np.sin(psi), np.cos(psi), 0], [0, 0, 1]]
    )
    rotation_y = np.array(
        [
            [np.cos(theta), 0, np.sin(theta)],
            [0, 1, 0],
            [-np.sin(theta), 0, np.cos(theta)],
        ]
    )
    rotation_x = np.array(
        [[1, 0, 0], [0, np.cos(phi), -np.sin(phi)], [0, np.sin(phi), np.cos(phi)]]
    )
    return rotation_x @ rotation_y @ rotation_z


def generate_dimer_geometry_euler_angles(geo, phi, theta, psi, translation):
    rotation_matrix = euler_angle_rotation(phi, theta, psi)
    second_monomer = []
    for atom in geo:
        pos = np.array(
            [float(atom.split()[1]), float(atom.split()[2]), float(atom.split()[3])]
        )
        new_pos = rotation_matrix @ pos
        new_pos[-1] += translation
        second_monomer.append(
            f"{atom.split()[0]} {new_pos[0]} {new_pos[1]} {new_pos[2]}"
        )
    out = [geo, second_monomer]
    if not check_no_overlaps(out):
        raise ValueError("Monomers overlap")
    return out
