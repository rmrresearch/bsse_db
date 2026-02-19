from pathlib import Path

ROOT = Path(__file__).resolve().parent

def make_loader_pml(xyz_dirname: str, pml_name: str):
    xyz_dir = ROOT / xyz_dirname
    xyz_files = sorted(xyz_dir.glob("*.xyz"), key=lambda p: int(p.stem.split("_", 1)[0]))

    
    if not xyz_files:
        raise FileNotFoundError(f"No .xyz files found in {xyz_dir}")

    pml_path = ROOT / pml_name

    with pml_path.open("w") as f:
        # Load each xyz file as its own object
        for xyz in xyz_files:
            obj_name = xyz.stem
            rel_path = xyz.relative_to(ROOT)  # so the pml works from ROOT
            f.write(f"load {rel_path}, {obj_name}\n")

        # Apply ball-and-stick to everything after loading
        f.write("preset.ball_and_stick(selection='all', mode=1)\n")

    print("WROTE:", pml_path)
    print("Objects:", len(xyz_files))

# Make translational and rotatinal .pml files (you can make translational too)
make_loader_pml("translational_pymol_xyz", "load_translational.pml")
make_loader_pml("rotational_pymol_xyz", "load_rotational.pml")
