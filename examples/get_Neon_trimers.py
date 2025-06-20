from bsse_calculator.generate_trimer_geometry import generate_trimer_geometry
from bsse_calculator.get_bsse_trimer import get_bsse_trimer
import json
import numpy as np
from tqdm import tqdm
import os
neon_monomer = [
    "Ne 0.00000000 0.00000000 0.00000000",
]
out = {}
for side_length in tqdm(list(np.arange(1.5, 4, 0.25))):
    try:
        trimer = generate_trimer_geometry(neon_monomer, side_length)
    except ValueError:
        print(f"Skipping {side_length}")
        continue
    os.makedirs(f"Ne_trimer_{side_length}", exist_ok=True)
    out[side_length] = get_bsse_trimer(trimer, f"Ne_trimer_{side_length}")

with open("data/Ne/Ne_trimer_energies_extra.json", "w") as f:
    json.dump(out, f)
