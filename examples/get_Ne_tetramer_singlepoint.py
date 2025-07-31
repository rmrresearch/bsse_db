from bsse_calculator.get_bsse_nmer import get_bsse_nmer
import json
import numpy as np
from tqdm import tqdm
import os

geometry = [
    ["Ne 0.00000000 0.00000000 -1.48101717"],
    ["Ne 1.75556579 1.75556579 0.49367239"],
    ["Ne 0.64258168 -2.39814747 0.49367239"],
    ["Ne -2.39814747 0.64258168  0.49367239"],
]
out = {}
os.makedirs(f"optimal_Ne_tetramer", exist_ok=True)

out["optimal"] = get_bsse_nmer(4, geometry, "optimal_Ne_tetramer", False)
with open("data/Ne/Ne_tetramer_singlepoint.json", "w") as f:
    json.dump(out, f)
