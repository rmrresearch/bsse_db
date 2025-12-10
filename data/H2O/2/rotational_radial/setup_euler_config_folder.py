#!/usr/bin/env python3
import argparse
import math
import os
import re
import shutil
import sys


def extract_angles_from_folder(folder_name: str):
    """
    From something like:
        'H20_dimer_0.0_0.0_1.5707963267948966_2.75'
    return the 4 numeric strings:
        ('0.0', '0.0', '1.5707963267948966', '2.75')

    Strategy:
    - split on '_'
    - find the token 'dimer'
    - everything after 'dimer' should be [alpha, beta, gamma, R_name]
    (R_name is the R encoded in the folder name; we keep it just for reference)
    """
    parts = folder_name.split("_")
    if "dimer" not in parts:
        raise ValueError(
            f"Folder name '{folder_name}' does not contain 'dimer' token."
        )

    idx = parts.index("dimer")
    suffix_parts = parts[idx + 1:]
    if len(suffix_parts) != 4:
        raise ValueError(
            f"Expected 4 numeric values after 'dimer' in '{folder_name}', "
            f"got {len(suffix_parts)}."
        )

    alpha_str, beta_str, gamma_str, r_name_str = suffix_parts
    return alpha_str, beta_str, gamma_str, r_name_str


def radians_str_to_degrees_str(rad_str: str) -> str:
    """
    Convert a string with radians to a nicely formatted degrees string.
    Example:
        '1.5707963267948966' -> '90.0'
        '0.0' -> '0.0'
    We use 1 decimal place to match your example.
    """
    val = float(rad_str)
    deg = math.degrees(val)
    return f"{deg:.1f}"


def compute_R_from_input_dimer(folder_path: str) -> float:
    """
    Read input_dimer.txt in folder_path, extract the two O atom coordinates
    from the 'geometry ... end' block, and compute the distance between them.

    Returns:
        R (float) = |r2 - r1|
    """
    input_path = os.path.join(folder_path, "input_dimer.txt")
    if not os.path.isfile(input_path):
        raise FileNotFoundError(
            f"input_dimer.txt not found in folder '{folder_path}'"
        )

    oxy_coords = []
    in_geom = False

    with open(input_path, "r") as f:
        for line in f:
            stripped = line.strip()
            lower = stripped.lower()

            if lower.startswith("geometry"):
                in_geom = True
                continue

            if in_geom:
                if lower == "end":
                    break

                if not stripped:
                    continue

                parts = stripped.split()
                # Expect lines like: O x y z
                if len(parts) >= 4 and parts[0] == "O":
                    try:
                        x = float(parts[1])
                        y = float(parts[2])
                        z = float(parts[3])
                        oxy_coords.append((x, y, z))
                    except ValueError:
                        continue  # skip malformed lines

    if len(oxy_coords) != 2:
        raise ValueError(
            f"Expected exactly 2 oxygen atoms in geometry block of '{input_path}', "
            f"found {len(oxy_coords)}."
        )

    (x1, y1, z1), (x2, y2, z2) = oxy_coords
    dx = x2 - x1
    dy = y2 - y1
    dz = z2 - z1

    R = math.sqrt(dx * dx + dy * dy + dz * dz)
    return R


def make_metadata_basename(source_basename: str,
                           alpha_deg: str,
                           beta_deg: str,
                           gamma_deg: str,
                           r_str: str) -> str:
    """
    Build the metadata base name using degrees and computed R:

    e.g. source_basename = 'H20_dimer_0.0_0.0_1.57_2.75'
         -> 'H20_dimer_0.0_0.0_90.0_2.76'
    """
    parts = source_basename.split("_")
    if "dimer" not in parts:
        raise ValueError(
            f"Folder name '{source_basename}' does not contain 'dimer' token."
        )
    idx = parts.index("dimer")
    # Keep everything up to and including 'dimer'
    prefix = "_".join(parts[:idx + 1])
    return f"{prefix}_{alpha_deg}_{beta_deg}_{gamma_deg}_{r_str}"


def confirm(prompt: str) -> bool:
    """Simple y/n confirmation."""
    ans = input(f"{prompt} (y/n): ").strip().lower()
    return ans == "y"


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Create numbered configuration folder for Euler-angle dimer "
            "and move/rename txt files. R is computed from input_dimer.txt."
        )
    )
    parser.add_argument("number", help="Identifier number (e.g. 1).")
    parser.add_argument(
        "source_folder",
        help="Source folder (e.g. H20_dimer_0.0_0.0_1.5707963267948966_2.75).",
    )
    args = parser.parse_args()

    number = str(args.number)
    source_folder = os.path.abspath(args.source_folder)

    if not os.path.isdir(source_folder):
        print(f"ERROR: Source folder does not exist: {source_folder}")
        sys.exit(1)

    source_basename = os.path.basename(os.path.normpath(source_folder))

    # ---- Extract angles & folder-encoded R (for reference) ----
    try:
        alpha_rad, beta_rad, gamma_rad, r_name_str = extract_angles_from_folder(
            source_basename
        )
    except ValueError as e:
        print(f"ERROR: {e}")
        sys.exit(1)

    # ---- Compute R from input_dimer.txt ----
    try:
        R_geom = compute_R_from_input_dimer(source_folder)
    except (FileNotFoundError, ValueError) as e:
        print(f"ERROR while computing R from input_dimer.txt: {e}")
        sys.exit(1)

    # Format computed R to 2 decimal places (string)
    R_str = f"{R_geom:.2f}"

    # Convert radians → degrees (as strings)
    alpha_deg = radians_str_to_degrees_str(alpha_rad)
    beta_deg = radians_str_to_degrees_str(beta_rad)
    gamma_deg = radians_str_to_degrees_str(gamma_rad)

    # New suffix: alpha_deg_beta_deg_gamma_deg_R (R from geometry)
    suffix = f"{alpha_deg}_{beta_deg}_{gamma_deg}_{R_str}"
    new_folder_name = f"{number}_{suffix}"
    new_folder_path = os.path.abspath(new_folder_name)

    txt_files = [f for f in os.listdir(source_folder) if f.lower().endswith(".txt")]

    # ----- INITIAL CHECK / SUMMARY -----
    print("\n=== SUMMARY ===")
    print(f"Source folder:              {source_folder}")
    print(f"New folder:                 {new_folder_path}")
    print(f"Parsed from name (radians): alpha={alpha_rad}, beta={beta_rad}, "
          f"gamma={gamma_rad}, R_name={r_name_str}")
    print(f"Converted (degrees):        alpha={alpha_deg}, beta={beta_deg}, "
          f"gamma={gamma_deg}")
    print(f"R from geometry (|r2-r1|):  R_geom={R_geom:.10f}  ->  R_str={R_str}")
    print(f"Text files found:           {len(txt_files)}")
    if txt_files:
        for f in txt_files:
            print(f"   - {f}")
    else:
        print("   (no .txt files found)")

    # Build the metadata filename using degrees + computed R
    metadata_basename = make_metadata_basename(
        source_basename, alpha_deg, beta_deg, gamma_deg, R_str
    )
    metadata_filename = f"{metadata_basename}.txt"

    print("\nThis will create:")
    print(f"  {new_folder_name}/CCSD_T/aug-cc-pvdz")
    print(f"  {new_folder_name}/CCSD_T/aug-cc-pvtz")
    print("and a metadata file:")
    print(f"  {metadata_filename}")
    print("containing:")
    print(
        f"  {number} --> alpha={alpha_deg}, beta={beta_deg}, "
        f"gamma={gamma_deg}, R={R_str}\n"
    )

    if not confirm("Proceed"):
        print("Aborted.")
        sys.exit(0)

    if os.path.exists(new_folder_path):
        print(f"ERROR: Target folder already exists: {new_folder_path}")
        sys.exit(1)

    # ----- CREATE FOLDER STRUCTURE -----
    ccsd_t_path = os.path.join(new_folder_path, "CCSD_T")
    aug_pvdz_path = os.path.join(ccsd_t_path, "aug-cc-pvdz")
    aug_pvtz_path = os.path.join(ccsd_t_path, "aug-cc-pvtz")

    os.makedirs(aug_pvdz_path)
    os.makedirs(aug_pvtz_path)

    # ----- METADATA FILE -----
    metadata_path = os.path.join(new_folder_path, metadata_filename)

    with open(metadata_path, "w") as f:
        f.write(
            f"{number} --> alpha={alpha_deg}, beta={beta_deg}, "
            f"gamma={gamma_deg}, R={R_str}\n"
        )

    print(f"Created metadata file: {metadata_path}")

    # ----- NEW AUTOMATIC RENAME RULES -----

    # Patterns for input_X_Y.txt and output_X_Y.txt
    patt_xy_input  = re.compile(r"input_(.+)_(.+)\.txt$")
    patt_xy_output = re.compile(r"output_(.+)_(.+)\.txt$")

    # Exact names for dimer/monomer
    dimer_in      = "input_dimer.txt"
    dimer_out     = "output_dimer.txt"
    monomer_in    = "input_monomer.txt"
    monomer_out   = "output_monomer.txt"

    moved_count = 0
    renamed_count = 0
    moved_details = []

    # ----- MOVE + AUTO-RENAME -----
    for fname in txt_files:
        src = os.path.join(source_folder, fname)
        new_name = None

        # 1. input_X_Y.txt → Y_X.nw
        m_in = patt_xy_input.match(fname)
        if m_in:
            X, Y = m_in.group(1), m_in.group(2)
            new_name = f"{Y}_{X}.nw"

        # 2. output_X_Y.txt → Y_X.out
        m_out = patt_xy_output.match(fname)
        if m_out:
            X, Y = m_out.group(1), m_out.group(2)
            new_name = f"{Y}_{X}.out"

        # 3–6. dimer/monomer fixed names
        if fname == dimer_in:
            new_name = "AB_AB.nw"
        elif fname == dimer_out:
            new_name = "AB_AB.out"
        elif fname == monomer_in:
            new_name = "A_A.nw"
        elif fname == monomer_out:
            new_name = "A_A.out"

        # If nothing matched, keep original name
        if new_name is None:
            new_name = fname
            print(f"Moving (no rename): {fname}")
        else:
            if new_name != fname:
                renamed_count += 1
                print(f"Moving & renaming: {fname}  ->  {new_name}")
            else:
                print(f"Moving (no rename): {fname}")

        dst = os.path.join(aug_pvdz_path, new_name)
        shutil.move(src, dst)
        moved_count += 1
        moved_details.append((fname, new_name))

    # ----- SUMMARY OF MOVED FILES -----
    print(f"\nMoved {moved_count} file(s) into {aug_pvdz_path}")
    print(f"Renamed {renamed_count} file(s) according to Euler renaming rules.\n")
    if moved_details:
        print("Files moved (original -> final name):")
        for old, new in moved_details:
            print(f"  {old}  ->  {new}")

    # ----- DELETE SOURCE FOLDER IF EMPTY -----
    remaining = os.listdir(source_folder)
    if len(remaining) == 0:
        print(f"\nSource folder '{source_folder}' is now empty.")
        if confirm("Delete this folder?"):
            os.rmdir(source_folder)
            print(f"Deleted folder: {source_folder}")
        else:
            print("Folder was NOT deleted.")
    else:
        print("\nSource folder still contains:")
        for item in remaining:
            print(f"  - {item}")
        print("Folder will NOT be deleted.")


if __name__ == "__main__":
    main()
