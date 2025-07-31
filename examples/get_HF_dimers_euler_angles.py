from bsse_calculator.get_geometry import generate_dimer_geometry_euler_angles
from bsse_calculator.get_bsse import get_bsses
import json
from tqdm import tqdm
from itertools import product
import os
import numpy as np
from collections import defaultdict

hf_monomer = ["H 0 0 -0.924", "F 0 0 0"]

out = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: defaultdict())))
for (phi, theta, psi), translation in tqdm(
    list(
        product(
            product(np.arange(0, 2 * np.pi, np.pi / 2), repeat=3),
            np.arange(1.75, 4, 0.5),
        )
    ),
):
    print(phi, theta, psi, translation)
    try:
        dimer = generate_dimer_geometry_euler_angles(
            hf_monomer, phi, theta, psi, translation
        )
    except ValueError:
        print(f"Skipping {phi}_{theta}_{psi}_{translation}")
        continue
    os.makedirs(f"HF_dimer_{phi}_{theta}_{psi}_{translation}", exist_ok=True)
    out[phi][theta][psi][translation] = get_bsses(
        dimer, f"HF_dimer_{phi}_{theta}_{psi}_{translation}", force_rerun=False
    )

with open("data/bsse_calculations_HF_dimer/HF_dimer_energies_euler.json", "w") as f:
    json.dump(out, f)
