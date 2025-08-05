from bsse_calculator.read_total_energy import read_energies
from pytest import approx
from test_bsse_helpers import get_assets_dir
import os

assets_dir = get_assets_dir()


def test_read_total_energy():
    output = read_energies(os.path.join(assets_dir, "sample_output.txt"))
    assert output["Total_energy"] == approx(-200.52055434739918, rel=1e-5)
    assert output["MP2_correlation_energy"] == approx(-0.459039477197523, rel=1e-5)
    assert output["CCSD_correlation_energy"] == approx(-0.0060180434527869475, rel=1e-5)
    assert output["SCF_energy"] == approx(-200.04640230404723, rel=1e-5)
