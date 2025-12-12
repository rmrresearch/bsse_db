#!/usr/bin/env python3
import os
import re
import argparse

# Patterns tailored to NWChem lines like:
# "         Total SCF energy =    -76.041603395490"
# " Total MP2 energy:            -76.266726270474308"
# " Total CCSD(T) energy:                 -76.279613385645391"

SCF_PATTERN = re.compile(
    r"Total\s+SCF\s+energy\s*[:=]\s*([\-0-9.+Ee]+)"
)
MP2_PATTERN = re.compile(
    r"Total\s+MP2\s+energy\s*[:=]\s*([\-0-9.+Ee]+)"
)
CCSDT_PATTERN = re.compile(
    r"Total\s+CCSD\(T\)\s+energy\s*[:=]\s*([\-0-9.+Ee]+)"
)


def extract_energy(lines, pattern):
    """
    Search through lines for the last occurrence of `pattern`.
    Returns a float or None.
    """
    value = None
    for line in lines:
        m = pattern.search(line)
        if m:
            try:
                value = float(m.group(1))
            except ValueError:
                continue
    return value


def parse_output_file(path):
    """Return (E_scf, E_mp2, E_ccsdt) from a single .out file."""
    try:
        with open(path, "r") as f:
            lines = f.readlines()
    except OSError as e:
        print(f"  [ERROR] Could not read {path}: {e}")
        return None, None, None

    E_scf = extract_energy(lines, SCF_PATTERN)
    E_mp2 = extract_energy(lines, MP2_PATTERN)
    E_ccsdt = extract_energy(lines, CCSDT_PATTERN)

    return E_scf, E_mp2, E_ccsdt


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Scan NWChem .out files and print SCF, MP2, CCSD(T) total "
            "energies and their correlation components."
        )
    )
    parser.add_argument(
        "root",
        nargs="?",
        default=".",
        help="Root directory to scan (default: current directory).",
    )
    args = parser.parse_args()

    root = os.path.abspath(args.root)
    print(f"Scanning for .out files under: {root}\n")

    any_found = False

    for dirpath, dirnames, filenames in os.walk(root):
        out_files = [f for f in filenames if f.lower().endswith(".out")]
        if not out_files:
            continue

        for fname in sorted(out_files):
            any_found = True
            full_path = os.path.join(dirpath, fname)
            rel_path = os.path.relpath(full_path, root)

            E_scf, E_mp2, E_ccsdt = parse_output_file(full_path)

            print(f"File: {rel_path}")

            # --- Totals ---
            if E_scf is not None:
                print(f"  SCF     total energy: {E_scf: .12f} hartree")
            else:
                print("  SCF     total energy: [not found]")

            if E_mp2 is not None:
                print(f"  MP2     total energy: {E_mp2: .12f} hartree")
            else:
                print("  MP2     total energy: [not found]")

            if E_ccsdt is not None:
                print(f"  CCSD(T) total energy: {E_ccsdt: .12f} hartree")
            else:
                print("  CCSD(T) total energy: [not found]")

            # --- Correlation pieces ---
            # MP2 correlation: E_mp2 - E_scf
            if E_scf is not None and E_mp2 is not None:
                Ecorr_mp2 = E_mp2 - E_scf
                print(f"  MP2 correlation energy:         {Ecorr_mp2: .12f} hartree")
            else:
                print("  MP2 correlation energy:         [cannot compute]")

            # CCSD(T) correlation *on top of* MP2: E_ccsdt - E_mp2
            if E_mp2 is not None and E_ccsdt is not None:
                Ecorr_ccsdt = E_ccsdt - E_mp2
                print(f"  CCSD(T) correlation component:  {Ecorr_ccsdt: .12f} hartree")
            else:
                print("  CCSD(T) correlation component:  [cannot compute]")

            print("-" * 60)

    if not any_found:
        print("No .out files found.")


if __name__ == "__main__":
    main()
