from bsse_calculator.get_geometry import generate_dimer_geometry_euler_angles
from bsse_calculator.get_bsse import get_bsses
import json
import numpy as np
from tqdm import tqdm
import os

hf_monomer = ["H 0 0 -0.924", "F 0 0 0"]
out = {}
for side_length in tqdm(list(np.arange(1.5, 4, 0.5))):
    try:
        dimer = generate_dimer_geometry_euler_angles(hf_monomer, 0, 0, 0, side_length)
    except ValueError:
        print(f"Skipping {side_length}")
        continue
    os.makedirs(f"HF_dimer_{side_length}_pvtz", exist_ok=True)
    out[side_length] = get_bsses(
        dimer, "aug-cc-pvtz", f"HF_dimer_{side_length}_pvtz", force_rerun=False
    )

with open("data/HF/HF_dimer_energies_pvtz.json", "w") as f:
    json.dump(out, f)
