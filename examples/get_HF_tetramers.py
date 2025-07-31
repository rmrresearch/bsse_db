from bsse_calculator.generate_nmer_geometry import generate_nmer_geometry
from bsse_calculator.get_bsse_nmer import get_bsse_nmer
import json
import numpy as np
from tqdm import tqdm
import os

hf_monomer = ["H 0 0 -0.924", "F 0 0 0"]
out = {}
for side_length in tqdm(list(np.arange(1.5, 4, 0.5))):
    tetramer = generate_nmer_geometry(hf_monomer, side_length, 4)
    os.makedirs(f"HF_tetramer_{side_length}", exist_ok=True)
    out[side_length] = get_bsse_nmer(
        4, tetramer, f"HF_tetramer_{side_length}", force=False
    )

with open("data/bsse_calculations_HF_dimer/HF_tetramer_energies.json", "w") as f:
    json.dump(out, f)
