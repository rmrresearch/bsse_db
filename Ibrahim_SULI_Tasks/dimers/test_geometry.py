from utils import geometry, dimer_input


def test_geometry():
  HF_bond_length = 1
  FF_distance = 4
  output = geometry(HF_bond_length, FF_distance)
  expected_output = [['H 0.0 0.0 -1', 'F 0.0 0.0 0.0'],
                     ['F 0.0 0.0 4', "H 0.0 0.0 5"]
  ]
  assert output == expected_output
