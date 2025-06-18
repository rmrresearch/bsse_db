from bsse_calculator.generate_dimer_geometry import generate_dimer_geometry
from bsse_calculator.get_bsse import get_bsses
import pandas as pd
import numpy as np
import os
import itertools
from tqdm import tqdm
import time

if __name__ == "__main__":
    outputs = []
    ammonia_monomer = [
        "N -0.00139164 -0.00902253 0.15230249",
        "H 0.84968456 -0.42148407 -0.23746780",
        "H -0.08064342 0.92116522 -0.26585262",
        "H -0.77845025 -0.55550557 -0.22700820",
    ]
    for z in tqdm(list(np.arange(1.5, 4, 0.25))):
        start_time = time.time()
        geometry = generate_dimer_geometry(ammonia_monomer, 0, 0, z)
        print(z)
        output = get_bsses(
            geometry,
            f"Ammonia_Ammonia_distance_{z}",
            force_rerun=False,
        )
        output["distance"] = z
        outputs.append(output)
        end_time = time.time()
        print(f"Time taken: {end_time - start_time}")
        print(output)
    os.makedirs("data", exist_ok=True)
    pd.DataFrame(outputs).to_csv(
        "data/BSSE_by_Ammonia_Ammonia_distance_nwchem.csv", index=False
    )
