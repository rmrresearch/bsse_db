from bsse_calculator.generate_dimer_geometry import generate_dimer_geometry
from bsse_calculator.get_bsse import get_bsse
import pandas as pd
import numpy as np
import time
water_monomer  = [
    "O -0.15218561 -0.00116398 0.00000000",
    "H 0.43776878 -0.76611046 0.00000000",
    "H 0.44759598 0.75607921 0.00000000",
]

outputs = []
for distance in np.arange(1.5, 6, .5):
    start = time.time()
    dimer_geometry = generate_dimer_geometry(water_monomer, distance)
    output = get_bsse(dimer_geometry, force_rerun=True)
    output["distance"] = distance
    outputs.append(output)
    print(f"Calculating BSSE for distance {distance} took {time.time() - start} seconds")

pd.DataFrame(outputs).to_csv("data/water_BSSE_by_water_water_distance_nwchem.csv", index=False)
