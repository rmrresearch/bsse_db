from bsse_calculator.get_bsse_nmer import get_bsse_nmer
import json
import numpy as np
from tqdm import tqdm
import os

geometry = [
    ["F 0 0 0.0", "H 0.0 0.947 0.0"],
    [
        "F 0.4187630618840785 2.4787902917832687 0.0",
        "H 1.3657630618840784 2.4787902917832687 0.0",
    ],
    [
        "F 2.8975533536673472 2.06002722989919 0.0",
        "H 2.8975533536673472 1.1130272298991903 0.0",
    ],
    [
        "F 2.4787902917832687 -0.4187630618840785 0.0",
        "H 1.5317902917832686 -0.4187630618840784 0.0",
    ],
]
out = {}
os.makedirs(f"optimal_HF_tetramer2", exist_ok=True)

out["optimal"] = get_bsse_nmer(4, geometry, "optimal_HF_tetramer2", False)
with open("data/bsse_calculations_HF_dimer/HF_tetramer_singlepoint2.json", "w") as f:
    json.dump(out, f)
