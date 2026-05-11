"""
make_pymol_radial_xyz_HF.py

Purpose
-------
Extracts the HF dimer geometry from each AB_AB.nw file in the radial
configurations and writes:

    1. pymol_radial_xyz/n_0.0_0.0_R.xyz  — one XYZ file per configuration
    2. load_radial.pml                    — PyMOL script to load all configurations

XYZ format:
    4
    HF dimer radial config n (0.0 0.0 R) aug-cc-pvtz
    H   0.000   0.000  -0.924
    F   0.000   0.000   0.000
    H   0.000   0.000   z1
    F   0.000   0.000   R

Only real atoms are written (H and F). Ghost functions (bqH, bqF) are
excluded.

Usage
-----
Run from HF/2/:

    python3 scripts/make_pymol_radial_xyz_HF.py

Output:
    HF/2/pymol_radial_xyz/    <- one .xyz per configuration
    HF/2/load_radial.pml      <- PyMOL loader script
"""

from pathlib import Path


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

ROOT      = Path(__file__).parent.parent   # HF/2/
RAD_DIR   = ROOT / "radial"
PYMOL_DIR = ROOT / "pymol_radial_xyz"
METHOD    = "CCSD_T"
BASIS     = "aug-cc-pvtz"


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def parse_geometry(nw_file: Path) -> list[tuple[str, str, str, str]]:
    """
    Parse real atom lines from the geometry block of a NWChem input file.
    Returns list of (element, x, y, z) tuples.
    Ghost functions (lines starting with 'bq') are excluded.
    """
    atoms = []
    in_geometry = False
    for line in nw_file.read_text().splitlines():
        stripped = line.strip()
        if stripped.lower() == "geometry":
            in_geometry = True
            continue
        if stripped.lower() == "end" and in_geometry:
            break
        if in_geometry and stripped:
            parts = stripped.split()
            if len(parts) == 4 and not parts[0].lower().startswith("bq"):
                atoms.append((parts[0], parts[1], parts[2], parts[3]))
    return atoms


def write_xyz(path: Path, config_name: str, atoms: list[tuple[str, str, str, str]]) -> None:
    """Write a standard XYZ file."""
    parts = config_name.split("_", 1)
    index = parts[0]
    coords = parts[1].replace("_", " ")
    lines = [
        str(len(atoms)),
        f"HF dimer radial config {index} ({coords}) aug-cc-pvtz",
    ]
    for elem, x, y, z in atoms:
        lines.append(f"{elem}  {x}  {y}  {z}")
    path.write_text("\n".join(lines) + "\n")


def write_pml(pml_path: Path, xyz_files: list[str]) -> None:
    """Write a PyMOL script that loads all XYZ files with ball-and-stick display."""
    lines = []
    for xyz_name in xyz_files:
        obj_name = Path(xyz_name).stem
        lines.append(f"load pymol_radial_xyz/{xyz_name}, {obj_name}")
    lines.append("hide everything")
    lines.append("show sticks")
    lines.append("show spheres")
    lines.append("set sphere_scale, 0.2")
    lines.append("set stick_radius, 0.3")
    lines.append("util.cbaw")
    lines.append("set antialias, 2")
    lines.append("set ray_opaque_background, off")
    lines.append("set stick_quality, 16")
    lines.append("set sphere_quality, 4")
    lines.append("set ray_trace_mode, 1")
    lines.append("zoom all")
    lines.append("orient")
    pml_path.write_text("\n".join(lines) + "\n")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:

    config_dirs = sorted([
        d for d in RAD_DIR.iterdir()
        if d.is_dir()
    ], key=lambda d: int(d.name.split("_")[0]))

    if not config_dirs:
        print(f"No configuration folders found in {RAD_DIR}")
        return

    PYMOL_DIR.mkdir(exist_ok=True)
    xyz_filenames = []
    ok = 0
    warn = 0

    for config_dir in config_dirs:
        config_name = config_dir.name
        nw_file = config_dir / METHOD / BASIS / "AB_AB.nw"

        if not nw_file.exists():
            print(f"  [WARN] Missing: {nw_file}")
            warn += 1
            continue

        atoms = parse_geometry(nw_file)
        if len(atoms) != 4:
            print(f"  [WARN] Expected 4 real atoms, got {len(atoms)} in {nw_file}")
            warn += 1
            continue

        xyz_name = f"{config_name}.xyz"
        xyz_path = PYMOL_DIR / xyz_name
        write_xyz(xyz_path, config_name, atoms)
        xyz_filenames.append(xyz_name)
        ok += 1

    # Write PyMOL loader script
    pml_path = ROOT / "load_radial.pml"
    write_pml(pml_path, sorted(xyz_filenames, key=lambda f: int(f.split("_")[0])))

    print(f"\nDone.")
    print(f"  {ok} XYZ files written to {PYMOL_DIR}")
    if warn:
        print(f"  {warn} warnings — check output above")
    print(f"  PyMOL script written to {pml_path}")
    print(f"\nTo visualize in PyMOL, run from HF/2/:")
    print(f"  pymol load_radial.pml")


if __name__ == "__main__":
    main()