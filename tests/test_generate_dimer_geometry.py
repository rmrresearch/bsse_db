from compare_geometry import compare_geometry
from bsse_calculator.generate_dimer_geometry import generate_dimer_geometry


def test_generate_dimer_geometry():
    input_geo = [
        "O -0.15218561 -0.00116398 0.00000000",
        "H 0.43776878 -0.76611046 0.00000000",
        "H 0.44759598 0.75607921 0.00000000",
    ]
    output_geo = [
        [
            "O -0.15218561 -0.00116398 0.00000000",
            "H 0.43776878 -0.76611046 0.00000000",
            "H 0.44759598 0.75607921 0.00000000",
        ],
        [
            "O 3.84781439 -0.00116398 0.0",
            "H 4.43776878 -0.76611046 0.0",
            "H 4.44759598 0.75607921 0.0",
        ],
    ]
    assert output_geo == generate_dimer_geometry(input_geo, 4)


def test_generate_dimer_geometry_with_multiple_translations():
    input_geo = [
        "O -0.15218561 -0.00116398 0.00000000",
        "H 0.43776878 -0.76611046 0.00000000",
        "H 0.44759598 0.75607921 0.00000000",
    ]
    output_geo = [
        [
            "O -0.15218561 -0.00116398 0.00000000",
            "H 0.43776878 -0.76611046 0.00000000",
            "H 0.44759598 0.75607921 0.00000000",
        ],
        [
            "O 1.84781439 0.99883602 3.0",
            "H 2.43776878 0.23388954 3.0",
            "H 2.44759598 1.75607921 3.0",
        ],
    ]
    compare_geometry(output_geo, generate_dimer_geometry(input_geo, 2, 1, 3))
