from bsse_calculator.generate_inputs_dimer import dimer_input
from test_helpers import get_assets_dir
import os

assets_dir = get_assets_dir()


def test_dimer_input(tmp_path):
    geometry = [["H 0 0 -1", "F 0 0 0"], ["F 0 0 4", "H 0 0 5"]]
    input_path = os.path.join(tmp_path, "input.txt")
    dimer_input(geometry, input_path)
    with (
            open(input_path, "r") as f,
            open(os.path.join(assets_dir, "input_dimer.txt"), "r") as
            f_expected,
    ):
        output = f.read()
        expected_output = f_expected.read()
        assert output == expected_output
