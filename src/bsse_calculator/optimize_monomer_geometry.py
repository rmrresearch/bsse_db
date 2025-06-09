from rdkit import Chem
from rdkit.Chem import AllChem
from bsse_calculator.prepare_input import prepare_input


def get_geometry(molecule_smile):
    mol = Chem.MolFromSmiles(molecule_smile)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, useRandomCoords=True)  # optimize with NwChem not rdkit
    conf = mol.GetConformer()
    geo = []
    for atom in mol.GetAtoms():
        print(atom.GetSymbol())
        pos = conf.GetAtomPosition(atom.GetIdx())
        geo.append(f"{atom.GetSymbol()} {pos.x:.2f} {pos.y:.2f} {pos.z:.2f}")
    return geo


def geometry_to_string(geometry):
    out = ""
    for line in geometry:
        out += "  " + line + "\n"
    return out


def optimize_monomer_geometry(molecule_smile, input_file_path):
    geometry = get_geometry(molecule_smile)
    prepare_input(geometry_to_string(geometry), input_file_path, "optimize")
