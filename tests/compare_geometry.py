from pytest import approx
def compare_geometry(geometry1, geometry2):
  for monomer1, monomer2 in zip(geometry1, geometry2):
    for line1, line2 in zip(monomer1, monomer2):
      symbol_1, x_1, y_1, z_1 = line1.split()
      symbol_2, x_2, y_2, z_2 = line2.split()
      assert symbol_1 == symbol_2
      assert float(x_1) == approx(float(x_2))
      assert float(y_1) == approx(float(y_2))
      assert float(z_1) == approx(float(z_2))