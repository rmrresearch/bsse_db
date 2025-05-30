import subprocess
import os
from rdkit import Chem
import py3Dmol
def get_geometry(HF_bond_length, FF_distance):
  """
  Returns geometry of dimer in format:
  [
    ['H 0.0 0.0 -1', 'F 0.0 0.0 0.0'],
    ['F 0.0 0.0 4', "H 0.0 0.0 5"]
  ]
  where the first list is monomer A and the second list is monomer B
  each string is of the form "element x y z"
  """
  return [
    [f'H 0.0 0.0 -{HF_bond_length}', 'F 0.0 0.0 0.0'],
    [f'F 0.0 0.0 {FF_distance}', f"H 0.0 0.0 {FF_distance + HF_bond_length}"]
  ]
def dimer_input(geometry, input_file_path):
  """
  Creates psi4 input file for geometry to calculate E(AB, AB)
  """
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
  """
  Creates psi4 input files for geometry to calculate E(A, AB) and E(B, AB)
  """
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
def read_total_energy(file_path):
  """
  Returns total energy from psi4 output file
  """
  with open(file_path, "r") as f:
    for line in f:
      if "Total Energy =" in line:
        return float(line.split()[-1])
    raise ValueError(f"Total energy not found in {file_path}")
def get_optimized_monomer_energy(geometry, input_file_path):
  """
  Creates psi4 input file for geometry to calculate E(A, A)
  """
  monomer = geometry[0]
  with open(input_file_path, "w") as f:
    f.write("memory 600 mb\n")
    f.write("\n")
    f.write("molecule HF {\n")
    for line in monomer:
      f.write(f"  {line}\n")
    f.write("}\n")
    f.write("\n")
    f.write("set basis aug-cc-pVDZ\n")
    f.write("energy('ccsd(t)')\n")
def get_bsse(geometry, scripts_dir =  "scripts") -> dict:
  """
  Calculates BSSE for geometry (of dimer)
  Returns a dictonary of the form
  {
    "Delta_E_AB_AB": float,
    "BSSE_A": float,
    "BSSE_B": float
  }
  where Delta_E_AB_AB = ΔE(AB, AB) and BSSE_A = ε(A, AB) and BSSE_B = ε(B, AB)
  """
  prev_dir = os.getcwd()
  os.makedirs(scripts_dir, exist_ok=True)
  os.chdir(scripts_dir)
  dimer_input(geometry, "input_dimer.txt")
  monomer_input(geometry, "input_A_AB.txt", "input_B_AB.txt")
  get_optimized_monomer_energy(geometry, "input_monomer.txt")
  subprocess.run(["psi4", "input_dimer.txt", "output_dimer.txt"])
  subprocess.run(["psi4", "input_A_AB.txt", "output_A_AB.txt"])
  subprocess.run(["psi4", "input_B_AB.txt", "output_B_AB.txt"])
  subprocess.run(["psi4", "input_monomer.txt", "output_monomer.txt"])
  total_energy = read_total_energy("output_dimer.txt")
  E_A_AB = read_total_energy("output_A_AB.txt")
  E_B_AB = read_total_energy("output_B_AB.txt")
  E_monomer = read_total_energy("output_monomer.txt")
  os.chdir(prev_dir)
  return {
    "Delta_E_AB_AB": total_energy - E_A_AB - E_B_AB,
    "BSSE_A": E_monomer - E_A_AB,
    "BSSE_B": E_monomer - E_B_AB
  }
def make_diagram(geometry, output_path):
  """
  Creates 3D diagram of geometry saved as html file
  """
  bonds = [
    (0, 1, Chem.BondType.SINGLE),
    (2, 3, Chem.BondType.SINGLE),
  ]
  mol = Chem.RWMol()
  for row in [*geometry[0], *geometry[1]]:
      atom_symbol = row.split()[0]
      mol.AddAtom(Chem.Atom(atom_symbol))
  for i, j, bond_type in bonds:
      mol.AddBond(i, j, bond_type)

  # Add conformer with coordinates
  conf = Chem.Conformer(mol.GetNumAtoms())
  for i, row in enumerate([*geometry[0], *geometry[1]]):
      _, x, y, z = row.split()
      conf.SetAtomPosition(i, Chem.rdGeometry.Point3D(float(x), float(y), float(z)))
  mol.AddConformer(conf)

  mol_block = Chem.MolToMolBlock(mol)

  view = py3Dmol.view(width=400, height=400)
  view.addModel(mol_block, 'mol')
  view.setStyle({'stick': {}, 'sphere': {'scale': 0.3}})
  view.zoomTo()
  html_str = view._make_html()
  with open(output_path, "w") as f:
    f.write(html_str)
