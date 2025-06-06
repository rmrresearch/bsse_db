import os

from bsse_calculator.generate_bsse_corrected_monomer_input import (
    bsse_corrected_monomer_input,
)

assets_dir = os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "test_assets", "nwchem"
)


def test_bsse_corrected_monomer_input(tmp_path):
    geometry = [["H 0 0 -1", "F 0 0 0"], ["F 0 0 4", "H 0 0 5"]]
    bsse_corrected_monomer_input(
        geometry,
        os.path.join(tmp_path, "input_A_AB.txt"),
        os.path.join(tmp_path, "input_B_AB.txt"),
    )
    for monomer in ["A", "B"]:
        with (
            open(os.path.join(tmp_path, f"input_{monomer}_AB.txt"), "r") as f,
            open(
                os.path.join(assets_dir, f"input_{monomer}_AB.txt"), "r"
            ) as f_expected,
        ):
            output = f.read()
            expected_output = f_expected.read()
            assert output == expected_output
