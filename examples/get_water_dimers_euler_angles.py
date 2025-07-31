from bsse_calculator.get_geometry import generate_dimer_geometry_euler_angles
from bsse_calculator.get_bsse import get_bsses
import json
from tqdm import tqdm
from itertools import product
import os
import numpy as np
from collections import defaultdict

water_monomer = [
    "O -0.15218561 -0.00116398 0.00000000",
    "H 0.43776878 -0.76611046 0.00000000",
    "H 0.44759598 0.75607921 0.00000000",
]
out = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: defaultdict())))
for (phi, theta, psi), translation in tqdm(
    list(product(np.arange(0.25, 4, 0.25), repeat=3)),
):
    try:
        dimer = generate_dimer_geometry_euler_angles(water_monomer, 0, 0, 0, ())
    except ValueError:
        print(f"Skipping {phi}_{theta}_{psi}_{translation}")
        continue
    os.makedirs(f"H20_dimer_{phi}_{theta}_{psi}_{translation}", exist_ok=True)
    out[phi][theta][psi][translation] = get_bsses(
        dimer, f"H20_dimer_{translation}", force_rerun=True
    )

with open("data/water/water_dimers_new.json", "w") as f:
    json.dump(out, f)
