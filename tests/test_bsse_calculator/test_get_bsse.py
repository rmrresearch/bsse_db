from bsse_calculator.get_bsse import get_bsses
from pytest import approx
from test_helpers import get_assets_dir
import os, shutil

assets_dir = get_assets_dir()


def test_get_bsse(tmp_path):
    geometry = [["H 0 0 -1", "F 0 0 0"], ["F 0 0 4", "H 0 0 5"]]
    for filename in [
            "output_dimer.txt",
            "output_A_AB.txt",
            "output_B_AB.txt",
            "output_monomer.txt",
    ]:
        src = os.path.join(assets_dir, filename)
        dst = os.path.join(tmp_path, filename)
        shutil.copyfile(src, dst)
    output = get_bsses(geometry, "aug-cc-pvdz", tmp_path, force_rerun=False)
    expected_output = {
        "Total_energy_Delta_E_AB_AB": 0.000515982992553,
        "Total_energy_BSSE_A": -3.5498e-11,
        "Total_energy_BSSE_B": -3.5498e-11,
    }
    for key in expected_output:
        assert output[key] == approx(expected_output[key], abs=1e-3)
