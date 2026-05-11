#!/usr/bin/env python3
import os

# Root folder where "rotational_radial" lives.
# Assuming you run this script from H2O/2
ROOT = "rotational_radial"

# Mapping of old -> new filenames
RENAME_MAP = {
    "AB_A.nw":  "A_AB.nw",
    "AB_A.out": "A_AB.out",
    "AB_B.nw":  "B_AB.nw",
    "AB_B.out": "B_AB.out",
}


def rename_in_dir(dirpath: str):
    """Rename AB_A/AB_B files inside a single aug-cc-pvdz directory."""
    files = os.listdir(dirpath)
    changed = False

    for old, new in RENAME_MAP.items():
        if old in files:
            old_path = os.path.join(dirpath, old)
            new_path = os.path.join(dirpath, new)

            if os.path.exists(new_path):
                print(f"  [WARN] {new} already exists in {dirpath}, skipping {old}")
                continue

            print(f"  Renaming: {old}  ->  {new}")
            os.rename(old_path, new_path)
            changed = True

    if not changed:
        print("  No matching AB_* files found here.")


def main():
    root_abs = os.path.abspath(ROOT)
    if not os.path.isdir(root_abs):
        print(f"ERROR: Root folder '{ROOT}' not found at {root_abs}")
        return

    print(f"Scanning under: {root_abs}\n")

    count_dirs = 0
    for dirpath, dirnames, filenames in os.walk(root_abs):
        # We only care about aug-cc-pvdz folders under CCSD_T
        base = os.path.basename(dirpath)
        if base == "aug-cc-pvdz" and "CCSD_T" in dirpath.split(os.sep):
            count_dirs += 1
            print(f"[{count_dirs}] Processing: {dirpath}")
            rename_in_dir(dirpath)
            print()

    if count_dirs == 0:
        print("No 'aug-cc-pvdz' directories under 'rotational_radial' were found.")
    else:
        print(f"Done. Processed {count_dirs} aug-cc-pvdz directory(ies).")


if __name__ == "__main__":
    main()
