"""
copy_rename_files_HF.py

Purpose
-------
Reorganize the HF dimer NWChem input/output files from the original
FF_distance_x_y_z/ flat structure into the standardized layout:

    translational/n_x_y_z/CCSD_T/aug-cc-pvdz/REAL_BASIS.nw
    translational/n_x_y_z/CCSD_T/aug-cc-pvdz/REAL_BASIS.out

where n is the configuration index (0-based, assigned by sorted
lexicographic order of the original folder names) and x, y, z are the
F-F displacement vector components in Angstroms.

Original files are preserved in HF_dimers/ as backup (the authoritative
backup is original_data/HF_dimers.tar).

File renaming convention (BASIS_REAL -> REAL_BASIS):
    input_E_A_A.txt   -> A_A.nw
    input_E_AB_AB.txt -> AB_AB.nw
    input_E_AB_A.txt  -> A_AB.nw
    input_E_AB_B.txt  -> B_AB.nw
    (same pattern for output_E_*.txt -> *.out)

Note on HF monomers:
    HF monomers A and B have the same geometry (both are HF molecules
    with H at z=-0.924 and F at z=0), just displaced in space.
    There is no separate B_B calculation — the parser uses A_A for
    both monomers when computing BSSE.

Usage
-----
Run from the HF/2/ directory:

    python3 scripts/copy_rename_files_HF.py

The script is safe to re-run: it skips destinations that already exist.
"""

import shutil
from pathlib import Path


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

ROOT      = Path(__file__).parent.parent   # HF/2/
SRC_DIR   = ROOT / "HF_dimers"
TRNL_DIR  = ROOT / "translational"
METHOD    = "CCSD_T"
BASIS     = "aug-cc-pvdz"
PREFIX    = "FF_distance_"


# ---------------------------------------------------------------------------
# File renaming map: original stem -> new stem
# (input_E_<key>.txt -> <value>.nw, output_E_<key>.txt -> <value>.out)
# ---------------------------------------------------------------------------

RENAME_MAP = {
    "A_A":   "A_A",
    "AB_AB": "AB_AB",
    "AB_A":  "A_AB",
    "AB_B":  "B_AB",
}


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def parse_coords(folder_name: str) -> tuple[str, str, str]:
    """Extract (x, y, z) strings from FF_distance_x_y_z folder name."""
    stem = folder_name.removeprefix(PREFIX)
    parts = stem.split("_")
    if len(parts) != 3:
        raise ValueError(f"Unexpected folder name format: {folder_name}")
    return parts[0], parts[1], parts[2]


def new_config_name(index: int, x: str, y: str, z: str) -> str:
    return f"{index}_{x}_{y}_{z}"


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:

    # 1. Verify source directory exists
    if not SRC_DIR.exists():
        print(f"ERROR: Source directory not found: {SRC_DIR}")
        print("Make sure HF_dimers/ exists in the HF/2/ directory.")
        return

    # 2. Collect and sort original configuration folders lexicographically
    orig_folders = sorted([
        f for f in SRC_DIR.iterdir()
        if f.is_dir() and f.name.startswith(PREFIX)
    ])

    if not orig_folders:
        print(f"No {PREFIX}* folders found in {SRC_DIR}")
        return

    print(f"Found {len(orig_folders)} configuration folders.")

    # 3. Create translational/ output directory
    TRNL_DIR.mkdir(exist_ok=True)

    # 4. Reorganize into translational/n_x_y_z/CCSD_T/aug-cc-pvdz/
    skipped = 0
    copied  = 0
    warned  = 0

    for index, folder in enumerate(orig_folders):
        try:
            x, y, z = parse_coords(folder.name)
        except ValueError as e:
            print(f"  [WARN] Skipping {folder.name}: {e}")
            warned += 1
            continue

        config_name = new_config_name(index, x, y, z)
        dest_dir = TRNL_DIR / config_name / METHOD / BASIS
        dest_dir.mkdir(parents=True, exist_ok=True)

        config_skipped = 0
        config_copied  = 0

        for orig_key, new_stem in RENAME_MAP.items():
            for ext, suffix in [("input", ".nw"), ("output", ".out")]:
                src = folder / f"{ext}_E_{orig_key}.txt"
                dst = dest_dir / f"{new_stem}{suffix}"

                if not src.exists():
                    print(f"  [WARN] Missing: {src.relative_to(ROOT)}")
                    warned += 1
                    continue

                if dst.exists():
                    config_skipped += 1
                    continue

                shutil.copy2(src, dst)
                config_copied += 1

        if config_skipped > 0 and config_copied == 0:
            print(f"  [SKIP] {config_name} — all files already exist")
            skipped += 1
        else:
            status = f"[{index:03d}]"
            if config_skipped > 0:
                status += f" (partial: {config_copied} copied, {config_skipped} skipped)"
            print(f"  {status} {folder.name} -> translational/{config_name}/")
            copied += 1

    # 5. Summary
    print(f"\n{'='*60}")
    print(f"Done.")
    print(f"  Configurations processed : {copied}")
    print(f"  Configurations skipped   : {skipped} (already existed)")
    print(f"  Warnings                 : {warned}")
    print(f"  Output directory         : {TRNL_DIR}")
    if warned == 0:
        print(f"\n✅ All files copied successfully.")
    else:
        print(f"\n⚠️  Some warnings were raised — check output above.")


if __name__ == "__main__":
    main()