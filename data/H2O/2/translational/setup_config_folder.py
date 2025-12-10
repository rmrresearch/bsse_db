#!/usr/bin/env python3
import argparse
import os
import re
import shutil
import sys


def extract_suffix_from_folder(folder_name: str) -> str:
    """
    From something like 'H2O_H2O_distance_0.5_1.5_1.75'
    return '0.5_1.5_1.75' (or last three tokens as fallback).
    """
    parts = folder_name.split("_")
    if "distance" in parts:
        idx = parts.index("distance")
        suffix_parts = parts[idx + 1:]
        if not suffix_parts:
            raise ValueError(
                f"Folder name '{folder_name}' has 'distance' but no values after it."
            )
        return "_".join(suffix_parts)

    if len(parts) < 3:
        raise ValueError(
            f"Folder name '{folder_name}' is too short to extract a suffix."
        )

    return "_".join(parts[-3:])


def parse_xyz_from_suffix(suffix: str):
    """
    Convert '0.5_1.5_1.75' → ('0.5', '1.5', '1.75').
    """
    parts = suffix.split("_")
    if len(parts) != 3:
        raise ValueError(f"Suffix '{suffix}' must contain exactly 3 numeric values.")
    return parts[0], parts[1], parts[2]


def confirm(prompt: str) -> bool:
    """Simple y/n confirmation."""
    ans = input(f"{prompt} (y/n): ").strip().lower()
    return ans == "y"


def main():
    parser = argparse.ArgumentParser(
        description="Create numbered configuration folder and move/rename txt files."
    )
    parser.add_argument("number", help="Identifier number (e.g. 41).")
    parser.add_argument("source_folder", help="Source folder containing txt files.")
    args = parser.parse_args()

    number = str(args.number)
    source_folder = os.path.abspath(args.source_folder)

    if not os.path.isdir(source_folder):
        print(f"ERROR: Source folder does not exist: {source_folder}")
        sys.exit(1)

    source_basename = os.path.basename(os.path.normpath(source_folder))

    try:
        suffix = extract_suffix_from_folder(source_basename)
    except ValueError as e:
        print(f"ERROR: {e}")
        sys.exit(1)

    new_folder_name = f"{number}_{suffix}"
    new_folder_path = os.path.abspath(new_folder_name)

    txt_files = [f for f in os.listdir(source_folder) if f.lower().endswith(".txt")]

    # ----- INITIAL CHECK / SUMMARY -----
    print("\n=== SUMMARY ===")
    print(f"Source folder:     {source_folder}")
    print(f"New folder:        {new_folder_path}")
    print(f"Suffix extracted:  {suffix}")
    print(f"Text files found:  {len(txt_files)}")
    if txt_files:
        for f in txt_files:
            print(f"   - {f}")
    else:
        print("   (no .txt files found)")

    print("\nThis will create:")
    print(f"  {new_folder_name}/CCSD_T/aug-cc-pvdz")
    print(f"  {new_folder_name}/CCSD_T/aug-cc-pvtz")
    print("and a metadata file:")
    print(f"  {source_basename}.txt")
    print("containing:")
    print(f"  {number} --> x=?, y=?, z=? extracted from suffix\n")

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
    x, y, z = parse_xyz_from_suffix(suffix)
    metadata_filename = f"{source_basename}.txt"
    metadata_path = os.path.join(new_folder_path, metadata_filename)

    with open(metadata_path, "w") as f:
        f.write(f"{number} --> x={x}, y={y}, z={z}\n")

    print(f"Created metadata file: {metadata_path}")

    # ----- AUTOMATIC RENAME PATTERNS -----
    input_pattern = re.compile(r"input_E_(.+_.+)\.txt")
    output_pattern = re.compile(r"output_E_(.+_.+)\.txt")

    moved_count = 0
    renamed_count = 0
    moved_details = []

    # ----- MOVE + AUTO-RENAME -----
    for fname in txt_files:
        src = os.path.join(source_folder, fname)

        new_name = fname  # default: keep name
        m_in = input_pattern.match(fname)
        m_out = output_pattern.match(fname)

        if m_in:
            xy = m_in.group(1)
            new_name = f"{xy}.nw"
        elif m_out:
            xy = m_out.group(1)
            new_name = f"{xy}.out"

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
    print(f"Renamed {renamed_count} file(s) based on input_E_/output_E_ patterns.\n")
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
