def geometry(HF_bond_length, FF_distance):
  return [
    [f'H 0.0 0.0 -{HF_bond_length}', 'F 0.0 0.0 0.0'],
    [f'F 0.0 0.0 {FF_distance}', f"H 0.0 0.0 {FF_distance + HF_bond_length}"]
  ]
def dimer_input(geometry, input_file_path):
  with open(input_file_path, "w") as f:
    f.write("memory 600 mb\n")
    f.write("\n")
    f.write("molecule HF {\n")
    for i, molecule in enumerate(geometry):
      for line in molecule:
        f.write(f"  {line}\n")
      if i != len(geometry) - 1:
        f.write("  --\n")
    f.write("}\n")
    f.write("\n")
    f.write("set basis aug-cc-pVDZ\n")
    f.write("energy('ccsd(t)')\n")
def monomer_input(geometry, input_A, input_B):
  with open(input_A, "w") as f:
    f.write("memory 600 mb\n")
    f.write("\n")
    f.write("molecule HF {\n")
    for i, molecule in enumerate(geometry):
      for line in molecule:
        f.write(f"  {line}\n" if i == 0 else f"  @{line}\n")
      if i != len(geometry) - 1:
        f.write("  --\n")
    f.write("}\n")
    f.write("\n")
    f.write("set basis aug-cc-pVDZ\n")
    f.write("energy('ccsd(t)')\n")
  with open(input_B, "w") as f:
    f.write("memory 600 mb\n")
    f.write("\n")
    f.write("molecule HF {\n")
    for i, molecule in enumerate(geometry):
      for line in molecule:
        f.write(f"  {line}\n" if i == 1 else f"  @{line}\n")
      if i != len(geometry) - 1:
        f.write("  --\n")
    f.write("}\n")
    f.write("\n")
    f.write("set basis aug-cc-pVDZ\n")
    f.write("energy('ccsd(t)')\n")
    