"""
copy_rename_files_HF_rot.py

Purpose
-------
Reorganize the HF dimer rotational NWChem input/output files from the
HF_files_extracted/ flat structure into the standardized layout:

    rotational/n_alpha_beta_gamma_R/CCSD_T/aug-cc-pvdz/REAL_BASIS.nw
    rotational/n_alpha_beta_gamma_R/CCSD_T/aug-cc-pvdz/REAL_BASIS.out

where:
    n           is the configuration index (0-based, assigned by
                lexicographic order of (alpha_deg, beta_deg, gamma_deg, R))
    alpha, beta, gamma  are the Euler angles in DEGREES
    R           is the F-F separation distance in Angstroms

The rotational dataset consists of 152 configurations recovered from
HF_files.tar (aug-cc-pVDZ basis set). The original folder names use
radians for the Euler angles:

    HF_dimer_alpha_beta_gamma_R
    e.g. HF_dimer_0.0_1.5707963267948966_0.0_2.75

These are converted to degrees in the config label:
    0.0                  -> 0.0
    1.5707963267948966   -> 90.0
    3.141592653589793    -> 180.0
    4.71238898038469     -> 270.0

Note on dataset completeness:
    The full rotational dataset has 304 configurations (4 alpha × 4 beta
    × 4 gamma × 5 R, minus 32 missing at alpha=0.0). Only 152 raw output
    files were recovered from the partially corrupted HF_files.tar —
    configs with alpha=180.0 and alpha=270.0 are missing. The complete
    energy data for all 304 configs is available in
    HF_dimer_energies_euler.json.

Monomer A geometry (fixed for all configs):
    H  0  0  -0.924
    F  0  0   0.000
Monomer B orientation is varied by the Euler angles.

Original files are preserved in HF_files_extracted/ as backup (the
authoritative backup is original_data/HF_files.tar).

File renaming convention:
    input_monomer.txt  -> A_A.nw       output_monomer.txt  -> A_A.out
    input_dimer.txt    -> AB_AB.nw     output_dimer.txt    -> AB_AB.out
    input_A_AB.txt     -> A_AB.nw      output_A_AB.txt     -> A_AB.out
    input_B_AB.txt     -> B_AB.nw      output_B_AB.txt     -> B_AB.out

Note on HF monomers:
    HF monomers A and B have the same bond geometry (H at z=-0.924,
    F at z=0). Monomer B is rotated by the Euler angles. There is no
    separate B_B calculation — the parser uses A_A for both monomers
    when computing BSSE.

Usage
-----
Run from the HF/2/ directory:

    python3 scripts/copy_rename_files_HF_rot.py

The script is safe to re-run: it skips destinations that already exist.
"""

import math
import shutil
from pathlib import Path


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

ROOT     = Path(__file__).parent.parent       # HF/2/
SRC_DIR  = ROOT / "HF_files_extracted"
ROT_DIR  = ROOT / "rotational"
METHOD   = "CCSD_T"
BASIS    = "aug-cc-pvdz"
PREFIX   = "HF_dimer_"


# ---------------------------------------------------------------------------
# File renaming map: original filename -> new filename
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

def rad_to_deg(rad_str: str) -> str:
    """
    Convert a radian string to a degrees string with 1 decimal place.
    e.g. '1.5707963267948966' -> '90.0'
         '0.0'                -> '0.0'
    """
    deg = math.degrees(float(rad_str))
    return f"{deg:.1f}"


def parse_folder(folder_name: str) -> tuple[str, str, str, str]:
    """
    Extract (alpha_rad, beta_rad, gamma_rad, R) from folder name.
    e.g. 'HF_dimer_0.0_1.5707963267948966_0.0_2.75'
      -> ('0.0', '1.5707963267948966', '0.0', '2.75')
    """
    stem = folder_name.removeprefix(PREFIX)
    # Split on '_' but we have floats with dots — need to be careful
    # Format is always: alpha_beta_gamma_R where each can be a float
    # Strategy: split on '_' and rejoin known float groups
    parts = stem.split("_")
    if len(parts) != 4:
        raise ValueError(
            f"Expected 4 numeric parts in '{folder_name}', got {len(parts)}: {parts}"
        )
    return parts[0], parts[1], parts[2], parts[3]


def sort_key(folder_name: str) -> tuple:
    """
    Sort key for lexicographic ordering by (alpha_deg, beta_deg, gamma_deg, R).
    """
    alpha, beta, gamma, r = parse_folder(folder_name)
    return (
        math.degrees(float(alpha)),
        math.degrees(float(beta)),
        math.degrees(float(gamma)),
        float(r),
    )


def new_config_name(index: int,
                    alpha_deg: str,
                    beta_deg: str,
                    gamma_deg: str,
                    r_str: str) -> str:
    return f"{index}_{alpha_deg}_{beta_deg}_{gamma_deg}_{r_str}"


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:

    # 1. Verify source directory exists
    if not SRC_DIR.exists():
        print(f"ERROR: Source directory not found: {SRC_DIR}")
        print("Make sure HF_files_extracted/ exists in the HF/2/ directory.")
        return

    # 2. Collect and sort original configuration folders
    orig_folders = sorted([
        f.name for f in SRC_DIR.iterdir()
        if f.is_dir() and f.name.startswith(PREFIX)
    ], key=sort_key)

    if not orig_folders:
        print(f"No {PREFIX}* folders found in {SRC_DIR}")
        return

    print(f"Found {len(orig_folders)} rotational configuration folders.")
    print()

    # 3. Create rotational/ output directory
    ROT_DIR.mkdir(exist_ok=True)

    # 4. Reorganize
    skipped = 0
    copied  = 0
    warned  = 0

    for index, folder_name in enumerate(orig_folders):
        try:
            alpha_rad, beta_rad, gamma_rad, r_str = parse_folder(folder_name)
        except ValueError as e:
            print(f"  [WARN] Skipping {folder_name}: {e}")
            warned += 1
            continue

        alpha_deg = rad_to_deg(alpha_rad)
        beta_deg  = rad_to_deg(beta_rad)
        gamma_deg = rad_to_deg(gamma_rad)

        config_name = new_config_name(index, alpha_deg, beta_deg, gamma_deg, r_str)
        dest_dir = ROT_DIR / config_name / METHOD / BASIS
        dest_dir.mkdir(parents=True, exist_ok=True)

        src_folder = SRC_DIR / folder_name
        config_skipped = 0
        config_copied  = 0

        for orig_name, new_name in RENAME_MAP.items():
            src = src_folder / orig_name
            dst = dest_dir / new_name

            if not src.exists():
                print(f"  [WARN] Missing: {folder_name}/{orig_name}")
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
                status += f" (partial: {config_copied} copied, "
                status += f"{config_skipped} skipped)"
            print(f"  {status} {folder_name} -> rotational/{config_name}/")
            copied += 1

    # 5. Summary
    print(f"\n{'='*60}")
    print(f"Done.")
    print(f"  Configurations processed : {copied}")
    print(f"  Configurations skipped   : {skipped} (already existed)")
    print(f"  Warnings                 : {warned}")
    print(f"  Output directory         : {ROT_DIR}")
    if warned == 0:
        print(f"\n✅ All files copied successfully.")
    else:
        print(f"\n⚠️  Some warnings were raised — check output above.")


if __name__ == "__main__":
    main()