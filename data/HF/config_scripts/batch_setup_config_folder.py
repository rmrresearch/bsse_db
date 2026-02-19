#!/usr/bin/env python3
import argparse
import glob
import os
import subprocess
import sys


def confirm(prompt: str) -> bool:
    ans = input(f"{prompt} (y/n): ").strip().lower()
    return ans == "y"


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Batch driver for setup_config_folder.py.\n"
            "It assigns consecutive numbers to many H2O_H2O_distance_* folders."
        )
    )
    parser.add_argument(
        "start_number",
        type=int,
        help="Starting integer for numbering (e.g. 76).",
    )
    parser.add_argument(
        "pattern",
        nargs="?",
        default="H2O_H2O_distance_*",
        help="Glob pattern for geometry folders (default: H2O_H2O_distance_*).",
    )
    parser.add_argument(
        "--script",
        default="setup_config_folder.py",
        help="Path to the single-folder script (default: setup_config_folder.py).",
    )

    args = parser.parse_args()

    start_number = args.start_number
    pattern = args.pattern
    script_path = args.script

    # Find matching folders
    candidates = glob.glob(pattern)
    folders = [f for f in candidates if os.path.isdir(f)]
    folders.sort()

    if not folders:
        print(f"No folders found matching pattern: {pattern}")
        sys.exit(1)

    print("The following folders will be processed:\n")
    current = start_number
    for folder in folders:
        print(f"  {current}  ->  {folder}")
        current += 1

    print()
    if not confirm("Proceed with these assignments"):
        print("Aborted.")
        sys.exit(0)

    # Run setup_config_folder.py for each folder
    current = start_number
    for folder in folders:
        print("\n" + "=" * 60)
        print(f"Processing folder: {folder}")
        print(f"Assigned number:   {current}")
        print("=" * 60)

        cmd = [sys.executable, script_path, str(current), folder]
        ret = subprocess.run(cmd)

        if ret.returncode != 0:
            print(f"\nWARNING: script returned non-zero exit code ({ret.returncode})")
            if not confirm("Continue with next folder?"):
                print("Stopping batch.")
                break

        current += 1

    print("\nBatch processing finished.")


if __name__ == "__main__":
    main()
