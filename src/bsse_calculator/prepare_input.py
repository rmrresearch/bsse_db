def prepare_input(geometry_str, input_file_path):
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
        f.write("task ccsd(t) energy\n")
