import msys
import numpy as np

import unittest
from openeye import oechem
from time import time
from contextlib import contextmanager
import tempfile

class Timer(object):
    def __init__(self):
        self.start()
    def start(self):
        self.bgn = time()
    def elapsed(self):
        return time() - self.bgn

@contextmanager
def TimedLogMessage(*args):
    t = Timer()
    yield
    #print("%.2f seconds to run:" % t.elapsed(), *args)

class Main(unittest.TestCase):

    def testFindDistinctBigFragments(self):
        mol = msys.Load('tests/files/1vcc.mae')
        mol.append(msys.Load('tests/files/2f4k.dms'))
        self.assertEqual(mol.ncts, 2)
        result = msys.FindDistinctFragments(mol, key='oechem_smiles')

    def testSmall(self):
        mol = msys.Load('tests/files/jandor-bad.sdf')

        # by default ConvertToOEChem doesn't assign the formal charges, should it?!
        contains_radicals = "C[C@H]([C@@H]1[C@H]2[C@H](C(=C(N2C1=O)C(=O)OCc3ccc(cc3)[N](=O)[O])Sc4ncccn4)OC)O"
        oemol = msys.ConvertToOEChem(mol)
        assert oemol.GetTitle() == '/u/nyc/gullingj/pub/examples/MMFF94_hypervalent/JANDOR.mae'
        self.assertEqual(oechem.OEMolToSmiles(oemol), contains_radicals)

        # Using OEChem to assign the formal charges
        properly_charged = "C[C@H]([C@@H]1[C@H]2[C@H](C(=C(N2C1=O)C(=O)OCc3ccc(cc3)[N+](=O)[O-])Sc4ncccn4)OC)O"
        oechem.OEAssignFormalCharges(oemol)
        self.assertEqual(oechem.OEMolToSmiles(oemol), properly_charged)

        # using msys to assign formal charges
        msys.AssignBondOrderAndFormalCharge(mol)
        amol = msys.AnnotatedSystem(mol)
        oemol = msys.ConvertToOEChem(mol)
        self.assertEqual(oechem.OEMolToSmiles(oemol), properly_charged)
        for matm, ratm in zip(mol.atoms, oemol.GetAtoms()):
            self.assertEqual(matm.atomic_number, ratm.GetAtomicNum())
            self.assertEqual(matm.formal_charge, ratm.GetFormalCharge())
            self.assertEqual(amol.valence(matm), ratm.GetExplicitValence())

    def testBig(self):
        mol = msys.Load('tests/files/2f4k.dms')
        msys.AssignBondOrderAndFormalCharge(mol)
        t = -time()
        oemol = msys.ConvertToOEChem(mol)
        assert oemol.GetDimension() == 3
        t += time()
        print("%s: %d atoms, %d bonds in %.3fs" % (mol.name, mol.natoms, mol.nbonds, t))

    def checkStereo(self, atom, expected):
        self.assertTrue(atom.IsChiral())
        self.assertTrue(atom.HasStereoSpecified())

        v = []
        for nbr in atom.GetAtoms():
            v.append(nbr)
        stereo = atom.GetStereo(v, oechem.OEAtomStereo_Tetrahedral)
        self.assertEqual(stereo, expected)

    def testIsotope(self):
        mol = msys.Load('tests/files/isotope.sdf')
        assert [a['isotope'] for a in mol.atoms] == [0,2,2]
        oemol = msys.ConvertToOEChem(mol)
        assert [a.GetIsotope() for a in oemol.GetAtoms()] == [0,2,2]
        mol = msys.ConvertFromOEChem(oemol)
        assert [a['isotope'] for a in mol.atoms] == [0,2,2]

    def testChiralAtoms(self):
        mol = msys.Load('tests/files/jandor.sdf')
        oemol = msys.ConvertToOEChem(mol)
        for r in oemol.GetAtoms():
            if r.GetIdx() in (0,4,7):
                self.checkStereo(r, oechem.OEAtomStereo_RightHanded)
            elif r.GetIdx() in (5,):
                self.checkStereo(r, oechem.OEAtomStereo_LeftHanded)

    def checkBondStereo(self, bond, expected):
        v = []
        for neigh in bond.GetBgn().GetAtoms():
            if neigh != bond.GetEnd():
                v.append(neigh)
                break
        for neigh in bond.GetEnd().GetAtoms():
            if neigh != bond.GetBgn():
                v.append(neigh)
                break
        stereo = bond.GetStereo(v, oechem.OEBondStereo_CisTrans)
        self.assertEqual(stereo, expected)

    def testBondStereo(self):
        sdf = 'tests/files/34106.sdf'
        mol = msys.Load(sdf)
        oemol = msys.ConvertToOEChem(mol)
        for bond in oemol.GetBonds():
            ai = bond.GetBgnIdx()
            aj = bond.GetEndIdx()
            if ai == 16 and aj == 17:
                self.assertTrue(bond.HasStereoSpecified(oechem.OEBondStereo_CisTrans))
                self.checkBondStereo(bond, oechem.OEBondStereo_Cis)
            else:
                self.assertFalse(bond.HasStereoSpecified(oechem.OEBondStereo_CisTrans))
                self.checkBondStereo(bond, oechem.OEBondStereo_Undefined)

    def testOmega(self):
        sdf = 'tests/files/methotrexate.sdf'
        mol = msys.Load(sdf)
        dmstmp = tempfile.NamedTemporaryFile(suffix='.dms')
        dms = dmstmp.name
        msys.Save(mol, dms)
        mol = msys.Load(dms)

        oemol = msys.ConvertToOEChem(mol)

        from openeye import oeomega
        opts = oeomega.OEOmegaOptions()
        opts.SetMaxConfs(1)
        opts.SetStrictStereo(False)
        omega = oeomega.OEOmega(opts)

        # OEOmegaOptions.SetCanonOrder canonicalizes the order of the
        # atoms in the molecule in order to make conformer generation
        # invariant of input atom order.
        for orig_idx, atom in enumerate(oemol.GetAtoms()):
            atom.SetData("orig_idx", orig_idx)
        assert omega(oemol)
        orig_idx_to_new_idx = {}
        for atom in oemol.GetAtoms():
            orig_idx_to_new_idx[atom.GetData("orig_idx")] = atom.GetIdx()

        for conf in oemol.GetConfs():
            coords = conf.GetCoords()
            npcrds = np.array([coords[orig_idx_to_new_idx[orig_idx]] for orig_idx in sorted(orig_idx_to_new_idx)])
            break

        cmol = mol.clone()
        cmol.positions = npcrds

        avg_length = 0.0
        for bond in cmol.bonds:
            vec = bond.first.pos - bond.second.pos
            avg_length += np.sqrt(vec.dot(vec))
        avg_length /= cmol.nbonds
        assert avg_length < 2.5, "average bond length is %.2f, atom order is likely screwed up!" % avg_length

        #msys.Save(cmol, 'meth_confs.sdf')

    def testConvertFromOEChem(self):
        testfiles = [
            'tests/files/jandor.sdf',
            'tests/files/test_UID_corrected.mol2',
            'tests/files/fused.sdf',
            'tests/files/3RYZ.pdb',
            'tests/files/1DUF.pdb',
            ]

        for testfile in testfiles:
            with oechem.oemolistream()(testfile) as ifs:
                mol = oechem.OEMol()
                with TimedLogMessage("OEReadMolecule"):
                    assert oechem.OEReadMolecule(ifs, mol)

                if 'pdb' in testfile:
                    oechem.OEAddExplicitHydrogens(mol)

                with TimedLogMessage("ConvertFromOEChem"):
                    system = msys.ConvertFromOEChem(mol)

                if 'pdb' in testfile:
                    assert len(system.residues) > 1
                    if '1DUF' in testfile:
                        assert len(system.chains) > 1

                with TimedLogMessage("ConvertToOEChem"):
                    new_oemol = msys.ConvertToOEChem(system)

                with TimedLogMessage("OEMolToSmiles"):
                    oechem.OEMolToSmiles(new_oemol)

                assert new_oemol.GetTitle() == mol.GetTitle()

                assert oechem.OEMolToSmiles(new_oemol) == oechem.OEMolToSmiles(mol)

    def testConvertAtoms(self):
        mol = msys.Load('tests/files/jandor.sdf')
        omol = msys.ConvertToOEChem(mol.atoms)
        assert omol.GetDimension() == 3
        new = msys.ConvertFromOEChem(omol)
        assert mol.positions.tolist() == new.positions.tolist()
        assert [a.atomic_number for a in mol.atoms] == [a.atomic_number for a in new.atoms]

        ifs = oechem.oemolistream()
        assert ifs.open('tests/files/jandor.sdf')
        oechem_mol = oechem.OEGraphMol()
        assert oechem.OEReadMolecule(ifs, oechem_mol)
        assert oechem.OEMolToSmiles(omol) == oechem.OEMolToSmiles(oechem_mol)

        # testing disconnected components
        atoms = mol.select('withinbonds 1 of hydrogen')
        omol = msys.ConvertToOEChem(atoms)
        new = msys.ConvertFromOEChem(omol)
        assert [a.pos.tolist() for a in atoms] == new.positions.tolist()
        assert [a.atomic_number for a in atoms] == [a.atomic_number for a in new.atoms]

    def testConvertFromOEChemFailures(self):
        mol = oechem.OEMol()

        # SMILES has no coordinates or hydrogens
        assert oechem.OESmilesToMol(mol, 'C1NCC[NH2+]C1')
        self.assertRaises(ValueError, msys.ConvertFromOEChem, mol)

        # This is a 2D SDF file
        with oechem.oemolistream()('tests/files/methotrexate.sdf') as ifs:
            assert oechem.OEReadMolecule(ifs, mol)
        self.assertRaises(ValueError, msys.ConvertFromOEChem, mol)

        # Suppressing hydrogens should cause ConvertFromOEChem to fail
        with oechem.oemolistream()('tests/files/fused.sdf') as ifs:
            assert oechem.OEReadMolecule(ifs, mol)
        assert oechem.OESuppressHydrogens(mol)
        self.assertRaises(ValueError, msys.ConvertFromOEChem, mol)

    def testNearlyPlanar(self):
        dms = msys.Load("tests/files/desresff420.dms")
        oe1 = msys.ConvertToOEChem(dms)
        oe2 = msys.ConvertToOEChem(dms.select("fragment 0"))
        s1 = oechem.OEMolToSmiles(oe1)
        s2 = oechem.OEMolToSmiles(oe2)
        assert s1 == s2

    def testLinearMoleculesWithNoHydrogens(self):
        from openeye import oeomega
        opts = oeomega.OEOmegaOptions()
        opts.SetMaxConfs(1)
        omega = oeomega.OEOmega(opts)

        for smi in ['O=S=O', 'c12c(onn1)onn2', 'C#N']:
            oemol = oechem.OEMol()
            assert oechem.OESmilesToMol(oemol, smi)
            assert omega(oemol)

            msys_sys = msys.ConvertFromOEChem(oemol)
            new_oemol = msys.ConvertToOEChem(msys_sys)

            self.assertEqual(new_oemol.GetDimension(), 3)
            self.assertEqual(oechem.OEMolToSmiles(new_oemol), smi)


if __name__=="__main__":
  unittest.main(verbosity=2)
