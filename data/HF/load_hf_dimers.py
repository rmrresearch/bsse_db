from glob import glob
import os
from pymol import cmd

# Folder that has all the .xyz files
XYZ_PATH = "/home/leothan/Documents/universities/ISU/PhD/Research/BSS_error_research/bsse_db/data/HF/HF_dimers_xyz"

def load_all_xyz():
    files = sorted(glob(os.path.join(XYZ_PATH, "*.xyz")))
    print("Found", len(files), "xyz files")

    if not files:
        print("No xyz files found in:", XYZ_PATH)
        return

    for f in files:
        obj_name = os.path.splitext(os.path.basename(f))[0]
        cmd.load(f, obj_name)
        cmd.hide("everything", obj_name)
        cmd.show("sticks", obj_name)

    print("Done. Each file is a separate object in the object panel.")

# Auto-run when you do "run load_h2o_dimers.py" in PyMOL
load_all_xyz()


