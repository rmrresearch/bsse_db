from bsse_calculator.get_bsse import get_bsse
from pytest import approx
import os, shutil

assets_dir = os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "test_assets", "nwchem"
)


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
    output = get_bsse(geometry, tmp_path, force_rerun=False)
    expected_output = {
        "Delta_E_AB_AB": 0.0010827957681129874,
        "BSSE_A": 0.00016099410720471496,
        "BSSE_B": 0.00016099410704839556,
    }
    for key in expected_output:
        assert output[key] == approx(expected_output[key], abs=1e-3)
