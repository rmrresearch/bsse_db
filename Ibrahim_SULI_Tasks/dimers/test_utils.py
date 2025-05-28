from utils import geometry, dimer_input
import os

def test_geometry():
  HF_bond_length = 1
  FF_distance = 4
  output = geometry(HF_bond_length, FF_distance)
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
def monomer_input(tmp_path):
  geometry = [
    ['H 0.0 0.0 -1', 
    'F 0.0 0.0 0.0'],
    ['F 0.0 0.0 4',
    "H 0.0 0.0 5"]
  ]
  dimer_input(geometry, os.path.join(tmp_path,"input_A.txt"), E_B_input = os.path.join(tmp_path,"input_B.txt"))
  for monomer in ["A", "B"]:
    with open(os.path.join(tmp_path,f"input_{monomer}.txt"), "r") as f, open(f"test_assets/input_{monomer}.txt", "r") as f_expected:
      output = f.read()
      expected_output = f_expected.read()
      assert output == expected_output