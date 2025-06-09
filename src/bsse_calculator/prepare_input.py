def prepare_input(geometry_str, input_file_path, task):
    """Takes a geometry string and writes it to a nwchem input file with name input_file_path.
    The task parameter is either "energy" or "optimize and determines whether to perform an single
    point energy calculation or geometry optimization"""
    with open(input_file_path, "w") as f:
        f.write("geometry\n")
        f.write(geometry_str)
        f.write("end\n")
        f.write("basis spherical\n")
        unique_atoms = set()
        for line in geometry_str.strip().split("\n"):
            atom = line.split()[0]
            if atom in unique_atoms:
                continue
            unique_atoms.add(atom)
            ghost = "bq" in atom
            if ghost:
                atom = atom.replace("bq", "")
                f.write(f"  bq{atom} library {atom} aug-cc-pvdz\n")
            else:
                f.write(f"  {atom} library aug-cc-pvdz\n")
        f.write("end\n")
        if task == "energy":
            f.write("task ccsd(t) energy\n")
        elif task == "optimize":
            f.write("task ccsd(t) optimize\n")
        else:
            raise ValueError(f"Unknown task: {task}")
