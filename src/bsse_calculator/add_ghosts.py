def add_ghosts(geo, real_atoms, basis_set):
    out = ""
    for monomer_idx, monomer in enumerate(geo):
        if monomer_idx in real_atoms:
            out += "\n".join(monomer)
            out += "\n"
        elif monomer_idx in basis_set:
            out += "\n".join([f"bq{line}" for line in monomer])
            out += "\n"
    return out
