from bsse_calculator.generate_dimer_geometry import generate_dimer_geometry


def test_generate_dimer_geometry():
    input_geo = [
        "O -0.15218561 -0.00116398 0.00000000",
        "H 0.43776878 -0.76611046 0.00000000",
        "H 0.44759598 0.75607921 0.00000000",
    ]
    generate_dimer_geometry(input_geo, 4)
    output_geo = [
        [
            "O -0.15218561 -0.00116398 0.00000000",
            "H 0.43776878 -0.76611046 0.00000000",
            "H 0.44759598 0.75607921 0.00000000",
        ],
        [
            "O -0.15218561 -0.00116398 4.0",
            "H 0.43776878 -0.76611046 4.0",
            "H 0.44759598 0.75607921 4.0",
        ],
    ]
    assert output_geo == generate_dimer_geometry(input_geo, 4)
