from bsse_calculator.generate_inputs_psi4 import GenerateInputsPsi4
import os, shutil
from pytest import approx

assets_dir = os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "test_assets", "psi4"
)
input_generator = GenerateInputsPsi4()


def test_geometry_FF_distance():
    HF_bond_length = 1
    FF_distance = 4
    output = input_generator.get_geometry(HF_bond_length, FF_distance=FF_distance)
    expected_output = [["H 0 0 -1", "F 0 0 0"], ["F 0 0 4", "H 0 0 5"]]
    assert output == expected_output


def test_geometry_coords():
    HF_bond_length = 1
    F2_coords = [1, 1, 1]
    output = input_generator.get_geometry(HF_bond_length, F2_coords=F2_coords)
    expected_output = [
        ["H 0 0 -1", "F 0 0 0"],
        ["F 1 1 1", "H 1 1 2"],
    ]
    assert output == expected_output


def test_dimer_input(tmp_path):
    geometry = [["H 0.0 0.0 -1", "F 0.0 0.0 0.0"], ["F 0.0 0.0 4", "H 0.0 0.0 5"]]
    input_path = os.path.join(tmp_path, "input.txt")
    input_generator.dimer_input(geometry, input_path)
    with (
        open(input_path, "r") as f,
        open(os.path.join(assets_dir, "input_dimer.txt"), "r") as f_expected,
    ):
        output = f.read()
        expected_output = f_expected.read()
        assert output == expected_output


def test_bsse_corrected_monomer_input(tmp_path):
    geometry = [["H 0.0 0.0 -1", "F 0.0 0.0 0.0"], ["F 0.0 0.0 4", "H 0.0 0.0 5"]]
    input_generator.bsse_corrected_monomer_input(
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


def test_monomer_input(tmp_path):
    geometry = [["H 0.0 0.0 -1", "F 0.0 0.0 0.0"], ["F 0.0 0.0 4", "H 0.0 0.0 5"]]
    input_generator.monomer_input(geometry, os.path.join(tmp_path, "input_monomer.txt"))
    with (
        open(os.path.join(tmp_path, f"input_monomer.txt"), "r") as f,
        open(os.path.join(assets_dir, "input_monomer.txt"), "r") as f_expected,
    ):
        output = f.read()
        expected_output = f_expected.read()
        assert output == expected_output


def test_read_total_energy():
    output = input_generator.read_energies(
        os.path.join(assets_dir, "sample_output.txt")
    )
    assert output["Total_energy"] == approx(-200.52055434739918, rel=1e-5)
    assert output["MP2_correlation_energy"] == approx(-0.4590394865732152, rel=1e-5)
    assert output["CCSD(T)_correlation_energy"] == approx(-0.0151125494, rel=1e-5)
    assert output["SCF_energy"] == approx(-200.04640230404723, rel=1e-5)


def test_get_bsse(tmp_path):
    geometry = [["H 0.0 0.0 -1", "F 0.0 0.0 0.0"], ["F 0.0 0.0 4", "H 0.0 0.0 5"]]
    for filename in [
        "output_dimer.txt",
        "output_A_AB.txt",
        "output_B_AB.txt",
        "output_monomer.txt",
    ]:
        src = os.path.join(assets_dir, filename)
        dst = os.path.join(tmp_path, filename)
        shutil.copyfile(src, dst)
    output = input_generator.get_bsse(geometry, tmp_path, force_rerun=False)
    expected_output = {
        "Delta_E_AB_AB": 0.0010827957681129874,
        "BSSE_A": 0.00016099410720471496,
        "BSSE_B": 0.00016099410704839556,
    }
    for key in expected_output:
        assert output[key] == approx(expected_output[key], abs=1e-3)
