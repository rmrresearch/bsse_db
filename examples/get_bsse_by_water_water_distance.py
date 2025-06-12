from bsse_calculator.generate_dimer_geometry import generate_dimer_geometry
from bsse_calculator.get_bsse import get_bsse
import pandas as pd
import numpy as np
import os
import itertools
from tqdm import tqdm

water_monomer = [
    "O -0.15218561 -0.00116398 0.00000000",
    "H 0.43776878 -0.76611046 0.00000000",
    "H 0.44759598 0.75607921 0.00000000",
]

if __name__ == "__main__":
    outputs = []
    water_monomer = [
        "O -0.15218561 -0.00116398 0.00000000",
        "H 0.43776878 -0.76611046 0.00000000",
        "H 0.44759598 0.75607921 0.00000000",
    ]
    for coords in tqdm(
        list(
            itertools.product(
                np.arange(0.25, 2, 0.25),
                np.arange(2, 3, 0.25),
                np.arange(0.25, 2, 0.25),
            )
        )
    ):
        if np.sqrt(coords[0] ** 2 + coords[1] ** 2 + coords[2] ** 2) < 2:
            continue
        geometry = generate_dimer_geometry(
            water_monomer, coords[0], coords[1], coords[2]
        )
        print(coords)
        output = get_bsse(
            geometry,
            f"H2O_H2O_distance_{'_'.join(map(str, coords))}",
            force_rerun=False,
        )
        output["x"] = coords[0]
        output["y"] = coords[1]
        output["z"] = coords[2]
        outputs.append(output)
    os.makedirs("data", exist_ok=True)
    pd.DataFrame(outputs).to_csv(
        "data/BSSE_by_H2O_H2O_distance_nwchem_more_distances.csv", index=False
    )
