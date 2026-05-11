"""
copy_rename_files_HF_pvtz.py

Purpose
-------
Reorganize the HF dimer radial NWChem input/output files from the
HF_pvtz_dimers.tar archive into the standardized layout:

    radial/n_0.0_0.0_R/CCSD_T/aug-cc-pvtz/REAL_BASIS.nw
    radial/n_0.0_0.0_R/CCSD_T/aug-cc-pvtz/REAL_BASIS.out

where n is the configuration index (0-based, assigned by ascending R
value) and R is the F-F separation distance in Angstroms along the
z-axis.

The radial dataset consists of 3 configurations extracted from
HF_pvtz_dimers.tar (aug-cc-pVTZ basis set):
    R = 2.5, 3.0, 3.5 Å

Monomer A geometry (fixed for all configs):
    H  0  0  -0.924
    F  0  0   0.000
Monomer B is displaced along the z-axis by R Angstroms.

Original tar is preserved in original_data/HF_pvtz_dimers.tar.

File renaming convention:
    input_monomer.txt  -> A_A.nw       output_monomer.txt  -> A_A.out
    input_dimer.txt    -> AB_AB.nw     output_dimer.txt    -> AB_AB.out
    input_A_AB.txt     -> A_AB.nw      output_A_AB.txt     -> A_AB.out
    input_B_AB.txt     -> B_AB.nw      output_B_AB.txt     -> B_AB.out

Note on HF monomers:
    HF monomers A and B have the same geometry (H at z=-0.924, F at
    z=0), just displaced along the z-axis. There is no separate B_B
    calculation — the parser uses A_A for both monomers when computing
    BSSE.

Note on config label format:
    The radial configs use the same n_x_y_z format as translational,
    with x=0.0, y=0.0, z=R. This is consistent with the Ne radial
    dataset convention and ensures parse_params() handles them
    correctly in build_raw_data.py.

Usage
-----
Run from the HF/2/ directory:

    python3 scripts/copy_rename_files_HF_pvtz.py

The script reads directly from the tar archive without extracting it.
It is safe to re-run: it skips destinations that already exist.
"""

import tarfile
from pathlib import Path, PurePosixPath


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

ROOT     = Path(__file__).parent.parent   # HF/2/
TAR_PATH = ROOT / "original_data" / "HF_pvtz_dimers.tar"
RAD_DIR  = ROOT / "radial"
METHOD   = "CCSD_T"
BASIS    = "aug-cc-pvtz"
PREFIX   = "HF_dimer_"
SUFFIX   = "_pvtz"


# ---------------------------------------------------------------------------
# File renaming map: original filename -> new stem + extension
# ---------------------------------------------------------------------------

RENAME_MAP = {
    "input_monomer.txt":  "A_A.nw",
    "output_monomer.txt": "A_A.out",
    "input_dimer.txt":    "AB_AB.nw",
    "output_dimer.txt":   "AB_AB.out",
    "input_A_AB.txt":     "A_AB.nw",
    "output_A_AB.txt":    "A_AB.out",
    "input_B_AB.txt":     "B_AB.nw",
    "output_B_AB.txt":    "B_AB.out",
}


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def parse_R(folder_name: str) -> str:
    """
    Extract R value string from folder name.
    e.g. 'HF_dimer_3.0_pvtz' -> '3.0'
    """
    stem = folder_name.removeprefix(PREFIX)
    r_str = stem.removesuffix(SUFFIX)
    return r_str


def new_config_name(index: int, r_str: str) -> str:
    """
    Build config label in n_x_y_z format with x=0.0, y=0.0, z=R.
    e.g. index=0, r_str='2.5' -> '0_0.0_0.0_2.5'
    """
    return f"{index}_0.0_0.0_{r_str}"


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:

    # 1. Verify tar exists
    if not TAR_PATH.exists():
        print(f"ERROR: Tar archive not found: {TAR_PATH}")
        return

    # 2. Discover config folders inside tar, sort by R value
    print(f"Reading tar: {TAR_PATH.name}")
    with tarfile.open(TAR_PATH, "r:") as tar:
        all_members = tar.getmembers()

    # Find unique top-level config directories
    config_folders = sorted(set(
        PurePosixPath(m.name).parts[0]
        for m in all_members
        if PurePosixPath(m.name).parts[0].startswith(PREFIX)
        and PurePosixPath(m.name).parts[0].endswith(SUFFIX)
    ), key=lambda f: float(parse_R(f)))

    if not config_folders:
        print(f"No {PREFIX}*{SUFFIX} folders found in tar.")
        return

    print(f"Found {len(config_folders)} radial configurations:")
    for f in config_folders:
        print(f"  {f}  (R = {parse_R(f)} Å)")
    print()

    # 3. Create radial/ output directory
    RAD_DIR.mkdir(exist_ok=True)

    # 4. Extract and reorganize
    skipped = 0
    copied  = 0
    warned  = 0

    with tarfile.open(TAR_PATH, "r:") as tar:
        for index, folder_name in enumerate(config_folders):
            r_str = parse_R(folder_name)
            config_name = new_config_name(index, r_str)
            dest_dir = RAD_DIR / config_name / METHOD / BASIS
            dest_dir.mkdir(parents=True, exist_ok=True)

            config_skipped = 0
            config_copied  = 0

            for orig_name, new_name in RENAME_MAP.items():
                # Build the internal tar path
                # Tar entries have leading './' prefix
                internal_path = f"./{folder_name}/{orig_name}"

                dst = dest_dir / new_name

                if dst.exists():
                    config_skipped += 1
                    continue

                # Extract file from tar
                try:
                    member = tar.getmember(internal_path)
                    f = tar.extractfile(member)
                    if f is None:
                        print(f"  [WARN] Could not read: {internal_path}")
                        warned += 1
                        continue
                    dst.write_bytes(f.read())
                    config_copied += 1
                except KeyError:
                    print(f"  [WARN] Missing in tar: {internal_path}")
                    warned += 1
                    continue

            if config_skipped > 0 and config_copied == 0:
                print(f"  [SKIP] {config_name} — all files already exist")
                skipped += 1
            else:
                status = f"[{index:02d}]"
                if config_skipped > 0:
                    status += f" (partial: {config_copied} copied, "
                    status += f"{config_skipped} skipped)"
                print(f"  {status} {folder_name} -> radial/{config_name}/")
                copied += 1

    # 5. Summary
    print(f"\n{'='*60}")
    print(f"Done.")
    print(f"  Configurations processed : {copied}")
    print(f"  Configurations skipped   : {skipped} (already existed)")
    print(f"  Warnings                 : {warned}")
    print(f"  Output directory         : {RAD_DIR}")
    if warned == 0:
        print(f"\n✅ All files copied successfully.")
    else:
        print(f"\n⚠️  Some warnings were raised — check output above.")


if __name__ == "__main__":
    main()