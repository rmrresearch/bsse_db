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
            "Batch driver for setup_euler_config_folder.py.\n"
            "Assigns consecutive numbers to many H20_dimer_* folders,\n"
            "processing them in batches of N (default 10)."
        )
    )
    parser.add_argument(
        "start_number",
        type=int,
        help="Starting integer for numbering (e.g. 1, 86, ...).",
    )
    parser.add_argument(
        "--pattern",
        default="H20_dimer_*",
        help="Glob pattern for Euler-angle folders (default: H20_dimer_*).",
    )
    parser.add_argument(
        "--script",
        default="setup_euler_config_folder.py",
        help="Path to the single-folder script (default: setup_euler_config_folder.py).",
    )
    parser.add_argument(
        "--batch-size",
        type=int,
        default=10,
        help="Number of folders to process per batch before asking to continue (default: 10).",
    )

    args = parser.parse_args()

    start_number = args.start_number
    pattern = args.pattern
    script_path = args.script
    batch_size = args.batch_size

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
        print(f"  {current:4d}  ->  {folder}")
        current += 1

    print()
    if not confirm("Proceed with these assignments"):
        print("Aborted.")
        sys.exit(0)

    # Process in batches
    current_number = start_number
    processed_in_batch = 0
    total = len(folders)

    for idx, folder in enumerate(folders):
        print("\n" + "=" * 70)
        print(f"Processing folder {idx+1}/{total}")
        print(f"  Assigned number: {current_number}")
        print(f"  Folder:          {folder}")
        print("=" * 70)

        cmd = [sys.executable, script_path, str(current_number), folder]
        ret = subprocess.run(cmd)

        if ret.returncode != 0:
            print(
                f"\nWARNING: {script_path} returned non-zero exit code ({ret.returncode}) "
                f"for folder '{folder}'."
            )
            if not confirm("Continue with next folder?"):
                print("Stopping batch run due to error.")
                break

        current_number += 1
        processed_in_batch += 1

        # After each batch_size folders, ask if we should continue
        if processed_in_batch >= batch_size and (idx + 1) < total:
            print(
                f"\nProcessed {processed_in_batch} folder(s) in this batch "
                f"(up to folder index {idx+1})."
            )
            if not confirm("Continue with the next batch?"):
                print("Stopping at this batch by user request.")
                break
            processed_in_batch = 0  # reset counter for the next batch

    print("\nBatch processing finished.")


if __name__ == "__main__":
    main()
