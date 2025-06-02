from utils import (
    get_geometry,
    dimer_input,
    monomer_input,
    read_total_energy,
    get_optimized_monomer_energy,
    get_bsse,
)
import os
from pytest import approx

assets_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "test_assets")

def test_get_geometry():
    HF_bond_length = 1
    FF_distance = 4
    output = get_geometry(HF_bond_length, FF_distance)
    expected_output = [
        ["H 0.0 0.0 -1", "F 0.0 0.0 0.0"],
        ["F 0.0 0.0 4", "H 0.0 0.0 5"],
    ]
    assert output == expected_output


def test_dimer_input(tmp_path):
    geometry = [["H 0.0 0.0 -1", "F 0.0 0.0 0.0"], ["F 0.0 0.0 4", "H 0.0 0.0 5"]]
    input_path = os.path.join(tmp_path, "input.txt")
    dimer_input(geometry, input_path)
    with open(input_path, "r") as f, open(
        os.path.join(assets_dir, "input_dimer.txt"), "r"
    ) as f_expected:
        output = f.read()
        expected_output = f_expected.read()
        assert output == expected_output


def test_monomer_input(tmp_path):
    geometry = [["H 0.0 0.0 -1", "F 0.0 0.0 0.0"], ["F 0.0 0.0 4", "H 0.0 0.0 5"]]
    monomer_input(
        geometry,
        os.path.join(tmp_path, "input_A_AB.txt"),
        os.path.join(tmp_path, "input_B_AB.txt"),
    )
    for monomer in ["A", "B"]:
        with open(os.path.join(tmp_path, f"input_{monomer}_AB.txt"), "r") as f, open(os.path.join(assets_dir, f"input_{monomer}_AB.txt"), "r"
        ) as f_expected:
            output = f.read()
            expected_output = f_expected.read()
            assert output == expected_output


def test_read_total_energy():
    output = read_total_energy(os.path.join(assets_dir, "sample_output.txt"))
    assert output == approx(-200.04640230404723, rel=1e-5)


def test_get_optimized_monomer_energy(tmp_path):
    geometry = [["H 0.0 0.0 -1", "F 0.0 0.0 0.0"], ["F 0.0 0.0 4", "H 0.0 0.0 5"]]
    get_optimized_monomer_energy(geometry, os.path.join(tmp_path, "input_monomer.txt"))
    with open(os.path.join(tmp_path, f"input_monomer.txt"), "r") as f, open(os.path.join(assets_dir, "input_monomer.txt"), "r"
    ) as f_expected:
        output = f.read()
        expected_output = f_expected.read()
        assert output == expected_output