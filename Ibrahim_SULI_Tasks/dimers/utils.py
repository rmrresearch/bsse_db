def geometry(HF_bond_length, FF_distance):
  return [
    [f'H 0.0 0.0 -{HF_bond_length}', 'F 0.0 0.0 0.0'],
    [f'F 0.0 0.0 {FF_distance}', f"H 0.0 0.0 {FF_distance + HF_bond_length}"]
  ]
