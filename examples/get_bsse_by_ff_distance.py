from bsse_calculator.generate_inputs_nwchem import get_geometry, get_bsse
import itertools
from tqdm import tqdm
import os
import pandas as pd
import numpy as np

if __name__ == "__main__":
    outputs = []
    for coords in tqdm(list(itertools.product(np.arange(0.25, 2, 0.25), repeat=3))):
        if np.sqrt(coords[0] ** 2 + coords[1] ** 2 + coords[2] ** 2) < 1.5:
            continue
        geometry = get_geometry(0.924, F2_coords=coords)
        print(coords)
        output = get_bsse(
            geometry, f"FF_distance_{'_'.join(map(str, coords))}", force_rerun=True
        )
        output["x"] = coords[0]
        output["y"] = coords[1]
        output["z"] = coords[2]
        print(output)
        outputs.append(output)
    os.makedirs("data", exist_ok=True)
    pd.DataFrame(outputs).to_csv("data/BSSE_by_FF_distance_nwchem.csv", index=False)
