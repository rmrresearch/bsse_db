from bsse_calculator.read_geometry_optimization_output import (
    read_geometry_optimization_output,
)
import os


assets_dir = os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "test_assets", "nwchem"
)


def test_read_geometry_optimization_output():
    out = read_geometry_optimization_output(
        os.path.join(assets_dir, "geometry_optimization_output.txt")
    )
    assert out == [
        "O -0.15218561 -0.00116398 0.00000000",
        "H 0.43776878 -0.76611046 0.00000000",
        "H 0.44759598 0.75607921 0.00000000",
    ]
