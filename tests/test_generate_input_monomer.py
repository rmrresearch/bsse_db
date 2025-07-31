from bsse_calculator.generate_inputs_monomer import monomer_input
import os

assets_dir = os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "test_assets", "nwchem"
)


def test_monomer_input(tmp_path):
    geometry = [["H 0 0 -1", "F 0 0 0"], ["F 0 0 4", "H 0 0 5"]]
    monomer_input(geometry, os.path.join(tmp_path, "input_monomer.txt"), "aug-cc-pvdz")
    with (
        open(os.path.join(tmp_path, f"input_monomer.txt"), "r") as f,
        open(os.path.join(assets_dir, "input_monomer.txt"), "r") as f_expected,
    ):
        output = f.read()
        expected_output = f_expected.read()
        assert output == expected_output
