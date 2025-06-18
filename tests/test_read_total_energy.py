from bsse_calculator.read_total_energy import read_energies
from pytest import approx
import os

assets_dir = os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "test_assets", "nwchem"
)


def test_read_total_energy():
    output = read_energies(os.path.join(assets_dir, "sample_output.txt"))
    assert output["Total_energy"] == approx(-200.52055434739918, rel=1e-5)
    assert output["MP2_correlation_energy"] == approx(-0.459039477197523, rel=1e-5)
    assert output["CCSD_correlation_energy"] == approx(-0.0060180434527869475, rel=1e-5)
    assert output["SCF_energy"] == approx(-200.04640230404723, rel=1e-5)
