import unittest

import sys
sys.path.append("../src")
from chemical_curation import curate

class TestMol(unittest.TestCase):

    def test_mol_from_smiles(self):
        smiles = "CCCCCC"
        mol = curate.Mol.from_smiles(smiles, precise_activities = {"kd": 13, "ic50": 34})
        print(mol)

    def test_mol_from_rdkit(self):
        from rdkit import Chem
        rdkit_mol = Chem.MolFromSmiles("CCCCCC")
        mol = curate.Mol.from_rdkit_mol(rdkit_mol, precise_activities = {"kd": 13, "ic50": 34})
        molblock = Chem.rdmolfiles.MolToMolBlock(rdkit_mol)
        print(molblock)
        print(mol)

    def test_mol_from_inchi(self):
        from rdkit import Chem
        inchi = "InChI=1S/C6H14/c1-3-5-6-4-2/h3-6H2,1-2H3"
        mol = curate.Mol.from_inchi(inchi, precise_activities = {"kd": 13, "ic50": 34})
        print(mol)

    def test_mol_from_mol_string(self):
        mol_string = '''
                     6  5  0  0  0  0  0  0  0  0999 V2000
                     0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
                     1.2990    0.7500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
                     2.5981   -0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
                     3.8971    0.7500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
                     5.1962   -0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
                     6.4952    0.7500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
                     1  2  1  0
                     2  3  1  0
                     3  4  1  0
                     4  5  1  0
                     5  6  1  0
                     M  END
                     '''
        mol = curate.Mol.from_mol_string(mol_string, precise_activities = {"kd": 13, "ic50": 34})
        print(mol)

    def test_mol_history(self):
        smiles = "CCCCCC"
        mol = curate.Mol.from_smiles(smiles, precise_activities = {"kd": 13, "ic50": 34})
        print(mol.history)

    def test_mixture(self):
        smiles = "CCCCCC.CC"
        mol = curate.Mol.from_smiles(smiles, precise_activities = {"kd": 13, "ic50": 34})
        self.assertTrue(mol.history.has_modification("Detected mixture"))

    #need to find example with salt that isn't detected as a mixture
    def test_salt(self):
        self.assertTrue(False)

    def test_organic(self):
        smiles = "OO"
        mol = curate.Mol.from_smiles(smiles, precise_activities = {"kd": 13, "ic50": 34})
        self.assertTrue(mol.history.has_modification("not organic"))
        self.assertTrue(mol.rejected)

    def test_illegal_atom(self):
        smiles = "C1CCC(CC1)[Mg]Br"
        mol = curate.Mol.from_smiles(smiles, precise_activities = {"kd": 13, "ic50": 34})
        self.assertTrue(mol.history.has_modification("allowed"))
        self.assertTrue(mol.rejected)


    #def test_metal(self):



if __name__ == '__main__':
    unittest.main()
