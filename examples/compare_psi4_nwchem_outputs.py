from bsse_calculator.generate_inputs_nwchem import GenerateInputsNwchem
from bsse_calculator.generate_inputs_psi4 import GenerateInputsPsi4

psi4_generator = GenerateInputsPsi4()
nwchem_generator = GenerateInputsNwchem()
geometry = psi4_generator.get_geometry(HF_bond_length=0.924, FF_distance=4.0)

data_dir = "/Users/ibrahims/Documents/Programming/undergrad_reasearch/ames_research_summer_2025/psi4/docstring_stuff/examples/data"


bsse_psi4 = psi4_generator.get_bsse(
    geometry,
    scripts_dir=f"{data_dir}/psi4",
    force_rerun=True,
)

bsse_nwchem = nwchem_generator.get_bsse(
    geometry,
    scripts_dir=f"{data_dir}/nwchem",
    force_rerun=True,
)

for key in bsse_psi4:
    print(
        f"{key} difference between psi4 and nwchem: {abs(bsse_psi4[key] - bsse_nwchem[key]) * 627.509}"
    )
