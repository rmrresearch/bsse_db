from bsse_calculator.generate_dimer_geometry import generate_dimer_geometry
from bsse_calculator.get_bsse import get_bsses
import json
import numpy as np
from tqdm import tqdm
import os

neon_monomer = [
    "Ne 0.00000000 0.00000000 0.00000000",
]
out = {}
for x, y, z in tqdm([[0.5, 0.5, 1.5]]):
    try:
        dimer = generate_dimer_geometry(neon_monomer, x, y, z)
    except ValueError:
        print(f"Skipping {x}, {y}, {z}")
        continue
    os.makedirs(f"Ne_dimer_{x}_{y}_{z}", exist_ok=True)
    out[(x, y, z)] = get_bsses(dimer, f"Ne_dimer_{x}_{y}_{z}")

# with open("data/Ne/Ne_dimer_energies_extra.json", "w") as f:
#     json.dump(out, f)
