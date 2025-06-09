import os
from bsse_calculator.optimize_monomer_geometry import optimize_monomer_geometry
import re

assets_dir = os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "test_assets", "nwchem"
)


def normalize_line(line):
    # Replace all numbers (including decimals and negatives) with a placeholder
    return re.sub(r"-?\d+\.?\d*", "<num>", line)


def test_optimize_monomer_geometry(tmp_path):
    optimize_monomer_geometry(
        "[OH2]",
        os.path.join(tmp_path, "input_optimize_water.txt"),
    )
    with (
        open(os.path.join(tmp_path, f"input_optimize_water.txt"), "r") as f,
        open(os.path.join(assets_dir, f"input_optimize_water.txt"), "r") as f_expected,
    ):
        for out_line, expected_line in zip(f, f_expected):
            out_line = normalize_line(out_line)
            expected_line = normalize_line(expected_line)
            assert (
                out_line == expected_line
            ), f"Lines differ:\n{out_line}\n{expected_line}"
