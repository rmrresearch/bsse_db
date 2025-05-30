from utils import get_geometry, dimer_input, monomer_input, read_energies, get_optimized_monomer_energy, get_bsse, make_diagram
import os, subprocess
from pytest import approx
import logging

def test_get_geometry():
  HF_bond_length = 1
  FF_distance = 4
  output = get_geometry(HF_bond_length, FF_distance)
  expected_output = [['H 0.0 0.0 -1', 'F 0.0 0.0 0.0'],
                     ['F 0.0 0.0 4', "H 0.0 0.0 5"]
  ]
  assert output == expected_output
def test_dimer_input(tmp_path):
  geometry = [
    ['H 0.0 0.0 -1', 
    'F 0.0 0.0 0.0'],
    ['F 0.0 0.0 4',
    "H 0.0 0.0 5"]
  ]
  input_path = os.path.join(tmp_path, "input.txt")
  dimer_input(geometry, input_path)
  with open(input_path, "r") as f, open("test_assets/input_dimer.txt", "r") as f_expected:
    output = f.read()
    expected_output = f_expected.read()
    assert output == expected_output
def test_monomer_input(tmp_path):
  geometry = [
    ['H 0.0 0.0 -1', 
    'F 0.0 0.0 0.0'],
    ['F 0.0 0.0 4',
    "H 0.0 0.0 5"]
  ]
  monomer_input(geometry, os.path.join(tmp_path,"input_A.txt"), os.path.join(tmp_path,"input_B.txt"))
  for monomer in ["A", "B"]:
    with open(os.path.join(tmp_path,f"input_{monomer}.txt"), "r") as f, open(f"test_assets/input_{monomer}_AB.txt", "r") as f_expected:
      output = f.read()
      expected_output = f_expected.read()
      assert output == expected_output
def test_read_energies():
  output = read_energies("test_assets/sample_output.txt")
  assert output["Total_energy"] == approx(-200.04640230404723, rel=1e-5)
  assert output["MP2_correlation_energy"] == approx(-0.4590394865732152, rel=1e-5)
  assert output["CCSD(T)_correlation_energy"] == approx(-0.0151125494 , rel=1e-5)
  assert output["SCF_energy"] == approx(-200.04640230404723, rel=1e-5)
def test_get_optimized_monomer_energy(tmp_path):
  geometry = [
    ['H 0.0 0.0 -1',
    'F 0.0 0.0 0.0'],
    ['F 0.0 0.0 4',
    "H 0.0 0.0 5"]
  ]
  get_optimized_monomer_energy(geometry, os.path.join(tmp_path,"input_monomer.txt"))
  with open(os.path.join(tmp_path,f"input_monomer.txt"), "r") as f, open(f"test_assets/input_monomer.txt", "r") as f_expected:
      output = f.read()
      expected_output = f_expected.read()
      assert output == expected_output
def test_calculate_bsse(tmp_path):
  geometry = [['H 0.0 0.0 -0.92', 'F 0.0 0.0 0.0'],
               ['F 0.0 0.0 2', 'H 0.0 0.0 2.92']]
  output = get_bsse(geometry, tmp_path)
  print(output)
  assert output["BSSE_A"] == approx(0.0004733027579248983, rel=1e-5)
  assert output["BSSE_B"] == approx(0.0004733027579248983, rel=1e-5)
  assert output["Delta_E_AB_AB"] == approx(0.023323579176050657, rel=1e-5)
def test_make_diagram(tmp_path):
  geometry = [['H 0.0 0.0 -0.92', 'F 0.0 0.0 0.0'],
               ['F 0.0 0.0 2', 'H 0.0 0.0 2.92']]
  make_diagram(geometry, os.path.join(tmp_path, "diagram.html"))
  logging.info(os.environ.get("OPEN_IMAGES") == "1")
  if os.environ.get("OPEN_IMAGES") == 1:
    subprocess.run(["open", os.path.join(tmp_path, "diagram.html")])