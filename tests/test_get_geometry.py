from bsse_calculator.get_geometry import get_geometry
import pytest


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
