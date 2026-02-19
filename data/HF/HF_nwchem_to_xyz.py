#!/usr/bin/env python3
import os

ROOT = "HF_dimers"         # where your many folders live
OUTPUT_DIR = "HF_dimers_xyz"
TARGET_NAME = "input_E_AB_AB.txt"


def extract_geometry(lines):
    """Extract atomic coordinates from NWChem geometry...end block."""
    geom = []
    in_geom = False

    for line in lines:
        low = line.strip().lower()

        if low.startswith("geometry"):
            in_geom = True
            continue

        if in_geom:
            if low == "end":
                break

            parts = line.split()
            if len(parts) >= 4:
                try:
                    elem = parts[0]
                    x, y, z = map(float, parts[1:4])
                    geom.append((elem, x, y, z))
                except ValueError:
                    pass

    return geom


def write_xyz(atoms, out_path, comment=""):
    """Write an XYZ file."""
    with open(out_path, "w") as f:
        f.write(f"{len(atoms)}\n")
        f.write(comment + "\n")
        for elem, x, y, z in atoms:
            f.write(f"{elem:2s} {x: .10f} {y: .10f} {z: .10f}\n")


def main():
    # Make sure output directory exists
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    for root, dirs, files in os.walk(ROOT):
        if TARGET_NAME in files:
            in_path = os.path.join(root, TARGET_NAME)
            folder_name = os.path.basename(root.rstrip(os.sep))

            # Output is ONE folder above
            out_path = os.path.join(OUTPUT_DIR, folder_name + ".xyz")

            with open(in_path, "r") as fh:
                lines = fh.readlines()

            atoms = extract_geometry(lines)

            if atoms:
                write_xyz(atoms, out_path, comment=in_path)
                print(f"[OK] {out_path}")
            else:
                print(f"[WARN] No geometry in {in_path}")


if __name__ == "__main__":
    main()
