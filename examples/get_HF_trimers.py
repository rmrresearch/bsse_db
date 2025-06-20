from bsse_calculator.generate_trimer_geometry import generate_trimer_geometry
from bsse_calculator.get_bsse_trimer import get_bsse_trimer
import json
import numpy as np
from tqdm import tqdm
import os
hf_monomer = ["H 0 0 -0.924", "F 0 0 0"]
out = {}
for side_length in tqdm(list(np.arange(1.5, 4, 0.25))):
    trimer = generate_trimer_geometry(hf_monomer, side_length)
    os.makedirs(f"HF_trimer_{side_length}", exist_ok=True)
    out[side_length] = get_bsse_trimer(trimer, f"HF_trimer_{side_length}")

with open("data/HF/HF_trimer_energies.json", "w") as f:
    json.dump(out, f)
