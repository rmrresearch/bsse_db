from bsse_calculator.utils import (
    read_energies,
    get_bsse,
    make_diagram,
)
import os, subprocess, shutil
from pytest import approx
from unittest.mock import patch

assets_dir = os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "test_assets", "psi4"
)


def test_read_total_energy():
    output = read_energies(os.path.join(assets_dir, "sample_output.txt"))
    assert output["Total_energy"] == approx(-200.04640230404723, rel=1e-5)
    assert output["MP2_correlation_energy"] == approx(-0.4590394865732152, rel=1e-5)
    assert output["CCSD(T)_correlation_energy"] == approx(-0.0151125494, rel=1e-5)
    assert output["SCF_energy"] == approx(-200.04640230404723, rel=1e-5)


def test_get_bsse(tmp_path):
    geometry = [["H 0.0 0.0 -1", "F 0.0 0.0 0.0"], ["F 0.0 0.0 4", "H 0.0 0.0 5"]]
    for filename in [
        "output_dimer.txt",
        "output_A_AB.txt",
        "output_B_AB.txt",
        "output_monomer.txt",
    ]:
        src = os.path.join(assets_dir, filename)
        dst = os.path.join(tmp_path, filename)
        shutil.copyfile(src, dst)
    output = get_bsse(geometry, tmp_path, force_rerun=False)
    expected_output = {
        "Delta_E_AB_AB": 0.0010827957681129874,
        "BSSE_A": 0.00016099410720471496,
        "BSSE_B": 0.00016099410704839556,
    }
    for key in expected_output:
        assert output[key] == approx(expected_output[key], rel=1e-5)


def test_make_diagram(tmp_path):
    geometry = [["H 0.0 0.0 -0.92", "F 0.0 0.0 0.0"], ["F 0.0 0.0 2", "H 0.0 0.0 2.92"]]
    make_diagram(geometry, os.path.join(tmp_path, "diagram.html"))
    if os.environ.get("OPEN_IMAGES"):
        subprocess.run(["open", os.path.join(tmp_path, "diagram.html")])
