from bsse_calculator.optimize_monomer_geometry import optimize_monomer_geometry
from bsse_calculator.read_geometry_optimization_output import (
    read_geometry_optimization_output,
)
import subprocess, os

if __name__ == "__main__":
    # os.makedirs("../data/ethene", exist_ok=True)
    # optimize_monomer_geometry("C=C", "../data/ethene/ethene.txt")
    # with open("../data/ethene/geometry_output.txt", "w") as f:
    #   subprocess.run(["nwchem", "../data/ethene/ethene.txt"], stdout=f)
    geo = read_geometry_optimization_output("data/ethene/geometry_output.txt")
    print(geo)
