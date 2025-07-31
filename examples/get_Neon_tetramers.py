from bsse_calculator.generate_nmer_geometry import generate_nmer_geometry
from bsse_calculator.get_bsse_nmer import get_bsse_nmer
import json
import numpy as np
from tqdm import tqdm
import os

Ne_monomer = [
    "Ne 0.00000000 0.00000000 0.00000000",
]
out = {}
for side_length in tqdm(list(np.arange(1.5, 4, 0.25))):
    try:
        tetramer = generate_nmer_geometry(Ne_monomer, side_length, 4)
    except ValueError:
        print(f"Skipping {side_length}")
        continue
    os.makedirs(f"Ne_tetramer_{side_length}", exist_ok=True)
    out[side_length] = get_bsse_nmer(
        4, tetramer, f"Ne_tetramer_{side_length}", force=False
    )

    os.makedirs("data/Ne", exist_ok=True)

with open("data/Ne/Ne_tetramer_energies.json", "w") as f:
    json.dump(out, f)
