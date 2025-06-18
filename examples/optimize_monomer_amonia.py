from bsse_calculator.optimize_monomer_geometry import optimize_monomer_geometry
from bsse_calculator.read_geometry_optimization_output import (
    read_geometry_optimization_output,
)

# optimize_monomer_geometry("N", "data/ammonia/nwchem_outputs/optimize_monomer.txt")

print(
    read_geometry_optimization_output(
        "data/ammonia/nwchem_outputs/optimize_monomer_output.txt"
    )
)
