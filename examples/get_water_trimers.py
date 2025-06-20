from bsse_calculator.generate_trimer_geometry import generate_trimer_geometry
from bsse_calculator.get_bsse_trimer import get_bsse_trimer
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
for side_length in tqdm(list(np.arange(1.5, 4, 0.25))):
    try:
        trimer = generate_trimer_geometry(water_monomer, side_length)
    except ValueError:
        print(f"Skipping {side_length}")
        continue
    os.makedirs(f"H20_trimer_{side_length}", exist_ok=True)
    out[side_length] = get_bsse_trimer(trimer, f"H20_trimer_{side_length}")

    os.makedirs("data/water", exist_ok=True)

with open("data/water/water_trimer_energies.json", "w") as f:
    json.dump(out, f)
