import openbabel as ob
from rdkit import Chem
from rdkit.Chem import AllChem
import rdkit
import pymol
file = open('data')
smiles = file.read().split('\n')

mol = Chem.MolFromSmiles(smiles[0])
mol = Chem.AddHs(mol)

mol_1 = Chem.MolFromSmiles(smiles[1])
mol_1 = Chem.AddHs(mol_1)




def minimize(mol):
    AllChem.EmbedMolecule(mol)
    prop = AllChem.MMFFGetMoleculeProperties(mol)
    ff = AllChem.MMFFGetMoleculeForceField(mol, prop)
    ff.Initialize()
    ff.Minimize(maxIts=1000000)
    print(ff.CalcEnergy())
    Chem.SanitizeMol(mol)
    return ff

ff = minimize(mol)
ff_1 = minimize(mol_1)


file = open('mol_pre.mol', 'w')
file.write(Chem.MolToMolBlock(mol, kekulize=False))
file.close()


file = open('mol_1_pre.mol', 'w')
file.write(Chem.MolToMolBlock(mol_1, kekulize=False))
file.close()

print(rdkit.Chem.rdMolAlign.AlignMol(mol, mol_1, atomMap=zip(range(mol.GetNumAtoms()), range(mol_1.GetNumAtoms()))))



file = open('mol.mol', 'w')
file.write(Chem.MolToMolBlock(mol, kekulize=True))
file.close()


file = open('mol_1.mol', 'w')
file.write(Chem.MolToMolBlock(mol_1, kekulize=True))
file.close()





