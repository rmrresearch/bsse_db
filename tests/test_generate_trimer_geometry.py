from bsse_calculator.generate_trimer_geometry import generate_trimer_geometry
import pytest
def test_generate_trimer_geometry():
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
            "O 1.8478143900000004 3.4629376351377545 0.00000000",
            "H 2.4377687800000003 2.6979911551377542 0.00000000",
            "H 2.4475959800000004 4.220180825137755 0.00000000",
        ],
        [
            "O 3.84781439 -0.00116398 0.00000000",
            "H 4.43776878 -0.76611046 0.00000000",
            "H 4.44759598 0.75607921 0.00000000",
        ],
    ]
    assert output_geo == generate_trimer_geometry(input_geo, 4)

def test_fails_with_overlap():
    input_geo = [
        "O -0.15218561 -0.00116398 0.00000000",
        "H 0.43776878 -0.76611046 0.00000000",
        "H 0.44759598 0.75607921 0.00000000",
    ]
    with pytest.raises(ValueError):
        generate_trimer_geometry(input_geo, 1.5)