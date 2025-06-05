from bsse_calculator.utils import get_bsse
from bsse_calculator.generate_inputs_psi4 import GenerateInputsPsi4
import itertools
from tqdm import tqdm
import os
import pandas as pd
import numpy as np

input_generator = GenerateInputsPsi4()
if __name__ == "__main__":
    bsse_data = pd.read_csv(os.path.join("data", "BSSE_by_FF_distance.csv"))
    files_needed = set(
        ["output_A_AB.txt", "output_B_AB.txt", "output_monomer.txt", "output_dimer.txt"]
    )
    for coords in tqdm(list(itertools.product(np.arange(0.25, 2, 0.25), repeat=3))):
        geometry = input_generator.get_geometry(0.924, F2_coords=coords)
        output = get_bsse(
            geometry, f"FF_distance_{'_'.join(map(str, coords))}", force_rerun=False
        )
        bsse_data = pd.concat(
            [
                bsse_data,
                pd.DataFrame(
                    {
                        "x": [coords[0]],
                        "y": [coords[1]],
                        "z": [coords[2]],
                        "BSSE_A": [output["BSSE_A"]],
                        "BSSE_B": [output["BSSE_B"]],
                    },
                ),
            ]
        )
    os.makedirs("data", exist_ok=True)
    bsse_data.to_csv(os.path.join("data", "BSSE_by_FF_distance.csv"), index=False)
