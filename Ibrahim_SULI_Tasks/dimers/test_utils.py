from utils import geometry, dimer_input


def test_geometry():
  HF_bond_length = 1
  FF_distance = 4
  output = geometry(HF_bond_length, FF_distance)
  expected_output = [['H 0.0 0.0 -1', 'F 0.0 0.0 0.0'],
                     ['F 0.0 0.0 4', "H 0.0 0.0 5"]
  ]
  assert output == expected_output
def test_dimer_input():
  geometry = [
    ['H 0.0 0.0 -1', 
    'F 0.0 0.0 0.0'],
    ['F 0.0 0.0 4',
    "H 0.0 0.0 5"]
  ]
  dimer_input(geometry, "input.txt")
  with open("input.txt", "r") as f, open("test_assets/input.txt", "r") as f_expected:
    output = f.read()
    expected_output = f_expected.read()
    assert output == expected_output