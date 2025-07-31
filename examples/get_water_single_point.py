from bsse_calculator.get_bsse_nmer import get_bsse_nmer_no_symmetry
import json
import numpy as np
from tqdm import tqdm
import os

geometry = [
    [
        "O 1.43578615 -1.81362614 0.02883928",
        "H 0.51948659 -2.05981895 0.27379913",
        "H 1.97323208 -1.99738588 0.80705304",
    ],
    [
        "O -0.20046914 0.37272269 -0.86419396",
        "H 0.63389482 -0.08242954 -0.65336297",
        "H -0.11854055 1.27366879 -0.50868257",
    ],
    [
        "O -1.33729724 -1.77663235 0.40232204",
        "H -1.94652043 -2.27621778 -0.15290166",
        "H -1.18633327 -0.92442897 -0.06636188",
    ],
    [
        "O 0.02509624 3.08759520 0.11490845",
        "H 0.17771767 3.32038150 1.03857842",
        "H 0.52221739 3.73343873 -0.40131992",
    ],
]
out = {}
os.makedirs(f"optimal_water_tetramer", exist_ok=True)
out["optimal"] = get_bsse_nmer_no_symmetry(4, geometry, "optimal_water_tetramer", False)
with open("data/water/water_tetramer_singlepoint.json", "w") as f:
    json.dump(out, f)
