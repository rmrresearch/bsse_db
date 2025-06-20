import numpy as np
from itertools import combinations, product
def generate_trimer_geometry(input_geo, side_length, force=False):
    x_trans = side_length * np.cos(np.pi / 3)
    y_trans = side_length * np.sin(np.pi / 3)
    monomer_1 = [
        f"{line.split()[0]} {float(line.split()[1]) + x_trans} {float(line.split()[2]) + y_trans} {line.split()[3]}"
        for line in input_geo
    ]
    # make equilateral triangle
    monomer_2 = [
        f"{line.split()[0]} {float(line.split()[1]) + side_length} {line.split()[2]} {line.split()[3]}"
        for line in input_geo
    ]
    out = [input_geo, monomer_1, monomer_2]
    if check_no_overlaps(out) or force:
        return [input_geo, monomer_1, monomer_2]
    else:
        raise ValueError("Overlapping atoms in trimer")

def check_no_overlaps(geometry):
    for monomers in combinations(geometry, 2):
        for atom1, atom2 in product(monomers[0], monomers[1]):
            pos1 = np.array([float(atom1.split()[1]), float(atom1.split()[2]), float(atom1.split()[3])])
            pos2 = np.array([float(atom2.split()[1]), float(atom2.split()[2]), float(atom2.split()[3])])
            if np.linalg.norm(pos1 - pos2) < 1.4:
                return False
    return True