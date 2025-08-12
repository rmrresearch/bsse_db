from bsse_calculator.get_geometry import (
    get_geometry,
    generate_dimer_geometry_euler_angles,
)
import pytest
import numpy as np


def test_geometry_FF_distance():
    HF_bond_length = 1
    FF_distance = 4
    output = get_geometry(HF_bond_length, FF_distance=FF_distance)
    expected_output = [["H 0 0 -1", "F 0 0 0"], ["F 0 0 4", "H 0 0 5"]]
    assert output == expected_output


def test_geometry_coords():
    HF_bond_length = 1
    F2_coords = [1, 1, 1]
    output = get_geometry(HF_bond_length, F2_coords=F2_coords)
    expected_output = [
        ["H 0 0 -1", "F 0 0 0"],
        ["F 1 1 1", "H 1 1 2"],
    ]
    assert output == expected_output


def test_bond_length_too_short():
    HF_bond_length = 0.5
    with pytest.raises(ValueError):
        get_geometry(HF_bond_length)


def test_coords_too_close():
    HF_bond_length = 1
    F2_coords = [0, 0.5, 0]
    with pytest.raises(ValueError):
        get_geometry(HF_bond_length, F2_coords=F2_coords)


def test_FF_distance_too_short():
    HF_bond_length = 1
    FF_distance = 1.3
    with pytest.raises(ValueError):
        get_geometry(HF_bond_length, FF_distance=FF_distance)


def test_geometry_euler_angles():
    water_monomer = [
        "O -0.15218561 -0.00116398 0.00000000",
        "H 0.43776878 -0.76611046 0.00000000",
        "H 0.44759598 0.75607921 0.00000000",
    ]
    expected_output = [
        [
            "O -0.15218561 -0.00116398 0.00000000",
            "H 0.43776878 -0.76611046 0.00000000",
            "H 0.44759598 0.75607921 0.00000000",
        ],
        [
            "O -0.15218561 -7.12732190635768e-20 3.99883602",
            "H 0.43776878 -4.691073613161532e-17 3.23388954",
            "H 0.44759598 4.6296499221417974e-17 4.75607921",
        ],
    ]
    output = generate_dimer_geometry_euler_angles(water_monomer, np.pi / 2, 0, 0, 4)
    assert output == expected_output
