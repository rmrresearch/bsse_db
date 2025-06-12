def generate_dimer_geometry(
    monomer_geometry, x_translation, y_translation=0, z_translation=0
):
    """
    Accepts a monomer and creates a dimer with the second atom translated by the parameters x, y, z
    """
    second_monomer = []
    for line in monomer_geometry:
        atom, x, y, z = line.split()
        second_monomer.append(
            f"{atom} {float(x)  + x_translation} {float(y) + y_translation} {float(z) + z_translation}"
        )
    return [monomer_geometry, second_monomer]
