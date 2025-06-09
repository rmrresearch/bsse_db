
def generate_dimer_geometry(monomer_geometry, distance):
  """
  Accepts a monomer and creates a dimer with the first atom in the geometry separated by distance
  """
  second_monomer = []
  for line in monomer_geometry:
    atom, x, y, z = line.split()
    second_monomer.append(f"{atom} {x} {y} {float(z) + distance}")
  return [monomer_geometry, second_monomer]