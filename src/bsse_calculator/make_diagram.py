import subprocess, os
from rdkit import Chem
import py3Dmol


def make_diagram(geometry, bonds, output_path):
    """
    Creates 3D diagram of geometry saved as html file
    """
    bonds = [(bond[0], bond[1], Chem.BondType.SINGLE) for bond in bonds]
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
    view.addModel(mol_block, "mol")
    view.setStyle({"stick": {}, "sphere": {"scale": 0.3}})
    view.zoomTo()
    html_str = view._make_html()
    with open(output_path, "w") as f:
        f.write(html_str)
