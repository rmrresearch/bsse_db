from bsse_calculator.get_geometry import generate_dimer_geometry_euler_angles
from bsse_calculator.get_bsse import get_bsses
import json
import numpy as np
from tqdm import tqdm
import os

water_monomer = [
    "O -0.15218561 -0.00116398 0.00000000",
    "H 0.43776878 -0.76611046 0.00000000",
    "H 0.44759598 0.75607921 0.00000000",
]
out = {}
for side_length in tqdm(list(np.arange(1.5, 4, 0.5))):
    try:
        dimer = generate_dimer_geometry_euler_angles(
            water_monomer, 0, 0, 0, side_length
        )
    except ValueError:
        print(f"Skipping {side_length}")
        continue
    os.makedirs(f"H20_dimer_{side_length}_pvtz", exist_ok=True)
    out[side_length] = get_bsses(
        dimer, "aug-cc-pvtz", f"H20_dimer_{side_length}_pvtz", force_rerun=False
    )

with open("data/H20/H20_dimer_energies_pvtz.json", "w") as f:
    json.dump(out, f)
