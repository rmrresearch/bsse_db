import numpy as np
from itertools import combinations, product


def regular_ngon_points(n, side_length=1.0, center=(0.0, 0.0)):
    """
    Generate the vertices of a regular n-gon.

    Parameters:
        n (int): Number of sides (vertices).
        radius (float): Radius of the circumcircle.
        center (tuple): (x, y) coordinates of the polygon's center.

    Returns:
        List of tuples: Vertices of the regular n-gon.
    """
    if n < 3:
        raise ValueError("n must be >= 3 for a valid polygon.")

    radius = (side_length / 2) / np.sin(np.pi / n)
    print(radius)
    cx, cy = center
    angles = np.linspace(0, 2 * np.pi, n, endpoint=False)
    return [
        (cx + radius * np.cos(theta), cy + radius * np.sin(theta)) for theta in angles
    ]


def generate_nmer_geometry(input_geo, side_length, n, force=False):
    ngon = regular_ngon_points(n, side_length)
    geo = [
        [
            f"{line.split()[0]} {float(line.split()[1]) + x_trans} {float(line.split()[2]) + y_trans} {line.split()[3]}"
            for line in input_geo
        ]
        for x_trans, y_trans in ngon
    ]
    if check_no_overlaps(geo) or force:
        print(geo)
        return geo
    else:
        raise ValueError("Overlapping atoms in trimer")


def check_no_overlaps(geometry):
    for monomers in combinations(geometry, 2):
        for atom1, atom2 in product(monomers[0], monomers[1]):
            pos1 = np.array(
                [
                    float(atom1.split()[1]),
                    float(atom1.split()[2]),
                    float(atom1.split()[3]),
                ]
            )
            pos2 = np.array(
                [
                    float(atom2.split()[1]),
                    float(atom2.split()[2]),
                    float(atom2.split()[3]),
                ]
            )
            if np.linalg.norm(pos1 - pos2) < 1.4:
                return False
    return True
