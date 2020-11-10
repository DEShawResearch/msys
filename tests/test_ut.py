import os, sys, unittest
import msys
from msys import knot, reorder, pfx, molfile, ff
from time import time
import numpy as NP
import json
import gzip
import shutil
import subprocess
import tempfile
import sqlite3
import random


def tmpfile(**kwds):
    return tempfile.NamedTemporaryFile(**kwds)


def vsize():
    cmd = "ps -p %d h -o vsz" % os.getpid()
    with os.popen(cmd) as p:
        s = p.read()
    return int(s)


def SaveJson(system, path, transform=None, maxDecimals=-1):
    data = msys.FormatJson(system, maxDecimals=maxDecimals).encode()
    if transform is not None:
        data = transform(data)
    with open(path, "wb") as ofs:
        ofs.write(data)


class TestJson(unittest.TestCase):
    def testExportJSONWithAllEmptyStringParameterTypeShouldNotProduceNamesField(self):
        s = msys.Load("tests/files/memo.dms")
        js = msys.FormatJson(s)
        assert "names" not in json.loads(js)
        new = msys.ParseJson(js)

    def testParsingJsonShouldNotCrash(self):
        msys.ParseJson(open("tests/files/crash.json").read())

    def test2f4k(self):
        tmp = tmpfile(suffix=".json")
        old = msys.Load("tests/files/2f4k.dms")
        msys.Save(old, tmp.name)
        new = msys.Load(tmp.name)
        self.assertEqual(old.hash(), new.hash())

    def testFormatParse(self):
        mol = msys.Load("tests/files/2f4k.dms")
        cell = mol.cell.flatten().tolist()
        js = msys.FormatJson(mol)
        d = json.loads(js)
        self.assertEqual(d["cell"], cell)
        new = msys.ParseJson(js)
        self.assertEqual(new.cell.flatten().tolist(), cell)

    def testCh4(self):
        old = msys.Load("tests/files/ch4.dms")
        js = msys.FormatJson(old)
        new = msys.ParseJson(js)
        ff = new.auxtable("viparr_info")
        assert ff.nparams == 1
        items = [(k, ff.param(0)[k]) for k in ff.props]
        self.assertEqual(
            items,
            [
                ("version", "3.0.3"),
                ("workdir", "/d/en/gullingj-0/p4/dms_inputs"),
                (
                    "cmdline",
                    "mae2dms"
                    " /proj/desres/root/Linux/x86_64/maeff_inputs/1.3.1/share/ch4.mae"
                    " ch4.dms ",
                ),
            ],
        )

    def testEmptyPosVelRes(self):
        mol = msys.CreateSystem()
        res = mol.addResidue()
        res.addAtom()
        res.addAtom()
        res.addAtom()
        mol.ct(0)["some"] = "data"
        js = msys.FormatJson(mol)
        print(len(js), js)
        d = json.loads(js)
        self.assertFalse("cell" in d)
        self.assertFalse("names" in d)
        self.assertFalse("chains" in d)
        self.assertFalse("residues" in d)
        self.assertFalse("pos" in d["i"])
        self.assertFalse("vel" in d["i"])
        self.assertFalse("res" in d["i"])
        self.assertFalse("name" in d["i"])
        self.assertFalse("tags" in d["i"])
        self.assertFalse("tags" in d["b"])
        self.assertTrue("c" in d)
        self.assertTrue("k" in d["c"][0])
        self.assertTrue("some" in d["c"][0]["k"])

        new = msys.ParseJson(js)
        assert new.ct(0)["some"] == "data"

    def testJsonIsReadable(self):
        mol = msys.ParseJson(open("tests/files/lig.json").read())

    def testJsonDeletedTermIsHandled(self):
        mol = msys.ParseJson(open("tests/files/lig.json").read())
        mol.table("exclusion").delTermsWithAtom(mol.atom(0))
        js = msys.FormatJson(mol)
        assert msys.ParseJson(js) is not None

    def testJsonSizes(self):
        import gzip
        import zlib
        import lzma
        import bz2

        compressors = [
            (None, lambda d: d, lambda d: d),
            (".gz", lambda d: gzip.compress(d, compresslevel=9), gzip.decompress),
            (".z", lambda d: zlib.compress(d, level=9), zlib.decompress),
            (".bz2", lambda d: bz2.compress(d, compresslevel=9), bz2.decompress),
        ]

        try:
            import zstandard

            def zstd(data):
                cctx = zstandard.ZstdCompressor(level=19)
                return cctx.compress(data)

            def unzstd(data):
                dctx = zstandard.ZstdDecompressor()
                return dctx.decompress(data)

            compressors.append((".zst", zstd, unzstd))
        except ImportError:
            pass

        compressors.append(
            (
                ".xz",
                lambda d: lzma.compress(d, format=lzma.FORMAT_XZ, preset=9),
                lzma.decompress,
            )
        )

        try:
            import brotli

            def br(data):
                return brotli.compress(data, quality=11)

            compressors.append((".br", br, brotli.decompress))
        except ImportError:
            pass

        sizes = []

        tmpdir = tempfile.TemporaryDirectory(prefix="msys_json")

        def test_all(name, mol, maxDecimals=-1):
            ref_terms = msys.Load("tests/files/cdk2-ligand-Amber14EHT.dms")
            ref_dms = tmpdir.name + "/" + name + ".dms"
            msys.SaveDMS(mol, ref_dms)
            these = [name]
            # print(name)
            for ext, transform, decompress in compressors:
                # print('\t', ext)
                fname = tmpdir.name + "/" + name + ".json"
                if ext is not None:
                    fname += ext
                SaveJson(mol, fname, transform=transform, maxDecimals=maxDecimals)
                these.append(os.path.getsize(fname))

                new = msys.ParseJson(decompress(open(fname, "rb").read()).decode())
                assert (
                    new.ct(0)["oesmi"]
                    == "C[NH2](C)C[C@@H](COc1ccc(cc1)Nc2cc(ncn2)Nc3c(cccc3F)F)O"
                )
                round_trip = tmpdir.name + "/" + name + "_from_" + str(ext) + ".dms"

                msys.SaveDMS(mol, round_trip)

                if "exclusion" not in new.table_names:
                    ff.experimental.build_pairs_and_exclusions(new)
                    places = 7
                    if maxDecimals != -1:
                        places = 0
                    assert_equal_nonbonded_terms(self, ref_terms, new, places=places)

                # these are slow and will ultimately fail once bond orders are removed since that relies on them
                # subprocess.check_call(['build/bin/dms-diff-ff', ref_dms, round_trip])
                # subprocess.check_call(['build/bin/dms-diff', ref_dms, round_trip])
            sizes.append(these)

        mol = msys.Load("tests/files/cdk2-ligand-Amber14EHT.dms")

        from openeye import oechem

        mol.ct(0)["oesmi"] = oechem.OEMolToSmiles(msys.ConvertToOEChem(mol))

        test_all("default", mol)
        mol.positions = NP.zeros(mol.positions.shape)
        test_all("no-pos", mol)
        mol.table("exclusion").remove()
        mol.table("pair_12_6_es").remove()
        test_all("no-excl", mol)

        def remove_bond_orders(mol):
            for atom in mol.atoms:
                atom.atomic_number = 0
                atom.formal_charge = 0
                atom.name = ""
            for bond in mol.bonds:
                bond.order = 0

        remove_bond_orders(mol)
        test_all("borders", mol)

        def remove_all_zero_term_properties(mol):
            for table in mol.tables:
                for name in table.term_props:
                    if all(t[name] == 0 for t in table.terms):
                        table.delTermProp(name)

        remove_all_zero_term_properties(mol)
        test_all("rmtprop", mol)

        test_all("maxdec", mol, maxDecimals=7)

        print()
        print("name", *(e for e, t, u in compressors), sep="\t")
        for these in sizes:
            print(*these, sep="\t")


def assert_equal_nonbonded_terms(self, mol, cpy, places=7):
    assert len(cpy.table("exclusion").terms) == len(mol.table("exclusion").terms)
    for term in cpy.table("exclusion").terms:
        atoms = [mol.atoms[a.id] for a in term.atoms]
        terms = mol.table("exclusion").findWithOnly(atoms)
        assert terms is not None
        assert len(terms) == 1
        assert [a.id for a in terms[0].atoms] == [a.id for a in term.atoms]

    keys = ["aij", "bij", "qij"]
    assert len(cpy.table("pair_12_6_es").terms) == len(mol.table("pair_12_6_es").terms)
    for term in cpy.table("pair_12_6_es").terms:
        atoms = [mol.atoms[a.id] for a in term.atoms]
        terms = mol.table("pair_12_6_es").findWithOnly(atoms)
        assert terms is not None
        assert len(terms) == 1
        assert [a.id for a in terms[0].atoms] == [a.id for a in term.atoms]

        for k in keys:
            self.assertAlmostEqual(term[k], terms[0][k], places=places)


class TestFF(unittest.TestCase):
    def test_build_pairs_and_exclusions(self):
        mol = msys.Load("tests/files/cdk2-ligand-Amber14EHT.dms")

        cpy = mol.clone()
        cpy.table("exclusion").remove()
        cpy.table("pair_12_6_es").remove()

        ff.experimental.build_pairs_and_exclusions(cpy)
        assert_equal_nonbonded_terms(self, mol, cpy)


class TestNeutralize(unittest.TestCase):
    def test1(self):
        from msys.neutralize import Neutralize

        mol = msys.Load("tests/files/2f4k.dms")
        msys.AssignBondOrderAndFormalCharge(mol)

        new = Neutralize(mol, verbose=True, concentration=0.1)
        print("old atoms", mol.natoms)
        print("new atoms:", new.natoms)
        self.assertEqual(len(new.selectIds("element Na")), 7)
        self.assertEqual(len(new.selectIds("element Cl")), 8)
        self.assertEqual(sum(a.formal_charge for a in new.atoms), 0)

        new = Neutralize(mol, verbose=True, concentration=0.2, cation="K")
        self.assertEqual(len(new.selectIds("element K")), 15)
        self.assertEqual(len(new.selectIds("element Cl")), 17)
        self.assertEqual(sum(a.formal_charge for a in new.atoms), 0)


class TestSmiles(unittest.TestCase):
    def testCharge(self):
        for sym, q in (
            ("", 0),
            ("+", 1),
            ("++", 2),
            ("-", -1),
            ("--", -2),
            ("+0", 0),
            ("+1", 1),
            ("+2", 2),
            ("-0", 0),
            ("-1", -1),
            ("-2", -2),
        ):
            smiles = "[P%s]" % sym
            mol = msys.FromSmilesString(smiles)
            self.assertEqual(mol.natoms, 1)
            fc = mol.atom(0).formal_charge
            self.assertEqual(q, fc, "%s: want %d got %d" % (smiles, q, fc))


class TestIndexedFile(unittest.TestCase):
    def testSdf(self):
        with tempfile.NamedTemporaryFile(suffix=".sdf") as tmp:
            shutil.copy("tests/files/cofactors.sdf", tmp.name)
            L = msys.IndexedFileLoader(tmp.name)
            self.assertEqual(L.path, tmp.name)
            self.assertEqual(len(L), 15)
            self.assertEqual(L[5].ct(0)["Name"], "NADP+")
            self.assertEqual(L[10].ct(0)["Name"], "FAD-CH2+")
            self.assertEqual(L[0].ct(0)["Name"], "dUMP anion")


class TestHash(unittest.TestCase):
    def testAtom(self):
        mol = msys.CreateSystem()
        h0 = mol.hash()
        mol.addAtom()
        h1 = mol.hash()
        mol.atom(0).charge = 0.1
        h2 = mol.hash()
        self.assertNotEqual(h0, h1)
        self.assertNotEqual(h0, h2)
        self.assertNotEqual(h1, h2)

        # add an atom, then remove it does not affect the hash
        # once we clone the system to remove dead residues, etc.
        mol.addAtom().remove()
        h3 = mol.clone().hash()
        self.assertEqual(h2, h3)

        m1 = mol.clone()
        m2 = mol.clone()
        m1.addAtomProp("foo", int)
        m2.addAtomProp("foo", float)
        self.assertNotEqual(m1.hash(), m2.hash())

    def testBond(self):
        mol = msys.CreateSystem()
        a1 = mol.addAtom()
        a2 = mol.addAtom()
        h0 = mol.hash()
        b12 = a1.addBond(a2)
        h1 = mol.hash()
        b12.order = 2
        h2 = mol.hash()
        mol.addBondProp("foo", str)
        h3 = mol.hash()
        self.assertEqual(len(set((h0, h1, h2, h3))), 4)

    def testTerm(self):
        mol = msys.Load("tests/files/small.mae")
        h = set()
        h.add(mol.hash())
        table = mol.table("stretch_harm")
        t = table.addTerm([mol.atom(0), mol.atom(2)], table.params.param(5))
        h.add(mol.hash())
        t["constrained"] = 42
        h.add(mol.hash())
        table.addTermProp("foo", int)
        h.add(mol.hash())
        table.params.addProp("bar", str)
        h.add(mol.hash())
        t["bar"] = "justinrocks"
        h.add(mol.hash())
        t["foo"] = 129
        h.add(mol.hash())
        self.assertEqual(len(h), 7)

    def testSystem(self):
        mol = msys.Load("tests/files/2f4k.dms")
        self.assertEqual(mol.hash(), 3583905226189957559)

    def testSorted(self):
        m1 = msys.CreateSystem()
        m1.name = "1"
        m1.addAtom()
        m1.addAtom()
        m1.updateFragids()
        tname = "stretch_harm"
        m1.addTableFromSchema(tname)

        m2 = m1.clone()
        m2.name = "2"

        # system 1: terms [0,1], [1,0], [0,1]
        table = m1.table(tname)
        param = table.params.addParam()
        atoms = m1.atoms
        table.addTerm(atoms, param)
        atoms.reverse()
        table.addTerm(atoms, param)
        atoms.reverse()
        table.addTerm(atoms, param)

        # system 1: terms [1,0], [0,1], [0,1]
        table = m2.table(tname)
        param = table.params.addParam()
        atoms = m2.atoms
        atoms.reverse()
        table.addTerm(atoms, param)
        atoms.reverse()
        table.addTerm(atoms, param)
        table.addTerm(atoms, param)

        self.assertEqual(m1.hash(sorted=True), m2.hash(sorted=True))
        self.assertNotEqual(m1.hash(sorted=False), m2.hash(sorted=False))

        # change parameters for first term in system 1 and third term in system 2
        # these systems are now equivalent, so we must break ties by parameters
        # or we'll get a different hash.
        m1.table(tname).term(0)["r0"] = 42
        m2.table(tname).term(2)["r0"] = 42
        self.assertEqual(m1.hash(sorted=True), m2.hash(sorted=True))
        self.assertNotEqual(m1.hash(sorted=False), m2.hash(sorted=False))


class TestHbond(unittest.TestCase):
    def test1(self):
        mol = msys.Load("tests/files/1vcc.mae", structure_only=True)
        finder = msys.HydrogenBondFinder(
            mol, "protein and name N", "protein and name O"
        )
        hbonds = finder.find(mol.positions)
        self.assertEqual(len(hbonds), 168)

        finder = msys.HydrogenBondFinder(
            mol, mol.select("protein and name N"), mol.select("protein and name O")
        )
        hbonds2 = finder.find(mol.positions)
        self.assertEqual(len(hbonds), 168)


class Contacts(unittest.TestCase):
    def compare(self, p1, p2):
        assert len(p1) > 0
        self.assertEqual(len(p1), len(p2))

        for e1, e2 in zip(p1, p2):
            self.assertEqual(e1[:2], e2[:2])
            d1 = e1[2]
            d2 = e2[2]
            self.assertTrue(abs(d1 - d2) < 1e-6)

    def testPairlist(self):
        mol = msys.Load("tests/files/1vcc.mae").clone("backbone and residue 0 to 5")
        pos = mol.positions.astype("f")
        pro = list(range(mol.natoms))
        atoms = [mol.atom(id) for id in pro]
        excl = msys.SpatialHash.Exclusions()
        eset = set()
        for tname in "stretch_harm", "angle_harm", "dihedral_trig":
            table = mol.table(tname)
            for t in table.findWithOnly(atoms):
                ids = [a.id for a in t.atoms]
                a, b = ids[0], ids[-1]
                excl.add(a, b)
                eset.add((a, b))
                eset.add((b, a))

        sh = msys.SpatialHash(pos, pro)
        i, j, d = sh.findContacts(3.0, pos, pro)
        old = sorted(
            p for p in zip(i, j, d) if p[0] < p[1] and tuple(p[:2]) not in eset
        )

        i, j, d = sh.findPairlist(3.0, excl)
        new = sorted(p for p in zip(i, j, d))
        self.assertEqual(old, new)

    def testNonbonded(self):
        """atom sets with no bonds in common"""
        mol = msys.Load("tests/files/2f4k.dms")
        pro = mol.selectIds("protein")
        wat = mol.selectIds("water")
        pos = mol.positions.astype("f")

        p1 = mol.findContactIds(3.2, pro, wat, pos)
        p1.sort()

        sh = msys.SpatialHash(pos, wat)
        i, j, d = sh.findContacts(3.2, pos, pro)
        p2 = list(zip(i, j, d))
        p2.sort()
        self.compare(p1, p2)

        i, j, d = sh.findContacts(3.2, pos, pro, reuse_voxels=True)
        p3 = list(zip(i, j, d))
        p3.sort()
        self.compare(p1, p3)

        sh.voxelize(8.0)
        i, j, d = sh.findContacts(3.2, pos, pro, reuse_voxels=True)
        p4 = list(zip(i, j, d))
        p4.sort()
        self.compare(p1, p4)

        sh.voxelize(1.0)
        i, j, d = sh.findContacts(3.2, pos, pro, reuse_voxels=True)
        p4 = list(zip(i, j, d))
        p4.sort()
        with self.assertRaises(AssertionError):
            self.compare(p1, p4)

        with self.assertRaises(RuntimeError):
            sh.voxelize(0)
        with self.assertRaises(RuntimeError):
            sh.voxelize(-1)

    def testSelf(self):
        """self-contacts"""
        mol = msys.Load("tests/files/jandor.sdf")
        pos = mol.positions.astype("f")
        p1 = mol.findContactIds(2.9)
        p1.sort()

        h = msys.SpatialHash(pos)
        i, j, d = h.findContacts(2.9, pos)
        p2 = zip(i, j, d)
        p2 = [
            (i, j, d)
            for i, j, d in p2
            if i < j and not mol.atom(i).findBond(mol.atom(j))
        ]
        p2.sort()
        self.compare(p1, p2)

    def testBorisov(self):
        mol = msys.Load("tests/files/spatial_hash_test.dms")
        id1 = 836
        id2 = 2264
        sel1 = "alpha and chain N"
        sel2 = "alpha and chain S"

        # check pbwithin selections
        pbsel = mol.selectIds("index %d and pbwithin 5 of index %d" % (id1, id2))
        self.assertEqual(pbsel, [id1])
        pbsel = mol.selectIds("index %d and pbwithin 5 of index %d" % (id2, id1))
        self.assertEqual(pbsel, [id2])

        pos = mol.positions.astype("f")
        box = mol.cell

        # check using spatial hash and no ids
        h = msys.SpatialHash([pos[id1]], box=box)
        found = h.findWithin(5, [pos[id2]])
        self.assertEqual(found, [0])

        h = msys.SpatialHash([pos[id2]], box=box)
        found = h.findWithin(5, [pos[id1]])
        self.assertEqual(found, [0])

        # check using spatial hash with single id
        h = msys.SpatialHash(pos, ids=[id1], box=box)
        found = h.findWithin(5, pos, ids=[id2])
        self.assertEqual(found, [id2])

        h = msys.SpatialHash(pos, ids=[id2], box=box)
        found = h.findWithin(5, pos, ids=[id1])
        self.assertEqual(found, [id1])

        # check using spatial hash with lots of ids
        idsel1 = mol.selectIds(sel1)
        idsel2 = mol.selectIds(sel2)
        self.assertTrue(id1 in idsel1)
        self.assertTrue(id2 in idsel2)

        h = msys.SpatialHash(pos, ids=idsel1, box=box)
        found = h.findWithin(5, pos, ids=idsel2)
        self.assertEqual(found, [id2])

        h = msys.SpatialHash(pos, ids=idsel2, box=box)
        found = h.findWithin(5, pos, ids=idsel1)
        self.assertEqual(found, [id1])

        # check findContacts
        h = msys.SpatialHash(pos, ids=idsel2, box=box)
        i1, i2, d = h.findContacts(5, pos, ids=idsel1)
        self.assertEqual(i1, [id1])
        self.assertEqual(i2, [id2])

    def testBadInputs(self):
        mol = msys.Load("tests/files/jandor.sdf")
        pos = mol.positions.astype("f")
        h = msys.SpatialHash(pos)
        with self.assertRaises(ValueError):
            h.findContacts(3, pos, [len(pos) + 40002])
        with self.assertRaises(ValueError):
            h.findWithin(3, pos, [len(pos) + 40002])

    def testSmallVoxels(self):
        pos = NP.zeros((5, 3), "f")
        pos[1] = [1000, 1000, 1000]
        pos[2] = [0.1, 0, 0]
        pos[3] = [1000, 1000, 1000.1]
        pos[4] = [100, -300, 100]
        h = msys.SpatialHash(pos[:2])
        r = h.findWithin(0.2, pos[2:]).tolist()
        self.assertEqual(r, [0, 1])

    def testPeriodicWater(self):
        mol = msys.Load("tests/files/2f4k.dms")
        pos = mol.positions.astype("f")
        os = mol.selectIds("oxygen and x>24")
        hs = mol.selectIds("hydrogen and x<-24")
        i1, j1, d1 = msys.SpatialHash(pos, hs, mol.cell).findContacts(2.5, pos, os)
        i2, j2, d2 = msys.SpatialHash(pos, os, mol.cell).findContacts(2.5, pos, hs)
        id1 = NP.stack((i1, j1), 1)
        id2 = NP.stack((i2, j2), 1)
        id3 = NP.stack((j2, i2), 1)
        assert len(id1) == 106
        self.assertEqual(id1.shape, id2.shape)
        m1 = dict(list(zip(list(map(tuple, id1)), d1)))
        m2 = dict(list(zip(list(map(tuple, id3)), d2)))
        self.assertEqual(sorted(m1.keys()), sorted(m2.keys()))
        for k, v1 in list(m1.items()):
            v2 = m2[k]
            self.assertTrue(abs(v1 - v2) < 1e-6)


class TestPropmap(unittest.TestCase):

    TABLE_NAME = "stretch_harm"

    def addprops(self, p):
        p["s"] = "32"
        p["i"] = -425
        p["f"] = 1.5

    def checkprops(self, p):
        self.assertEqual(type(p["s"]), str)
        self.assertEqual(type(p["i"]), int)
        self.assertEqual(type(p["f"]), float)

    def setUp(self):
        self.m = msys.CreateSystem()
        self.table = self.m.addTableFromSchema(self.TABLE_NAME)

    def testKeys(self):
        p = self.table.props
        self.assertEqual(list(p.keys()), [])
        p["foo"] = 1
        self.assertEqual(list(p.keys()), ["foo"])
        del p["foo"]
        self.assertEqual(list(p.keys()), [])

    def testTypes(self):
        p = self.table.props
        self.addprops(p)
        self.checkprops(p)

    def testChangeType(self):
        p = self.table.props
        p["foo"] = 32
        self.assertEqual(p["foo"], 32)
        p["foo"] = 32.1
        self.assertEqual(p["foo"], 32.1)
        p["foo"] = "32.1"
        self.assertEqual(p["foo"], "32.1")
        p["foo"] = 32
        self.assertEqual(p["foo"], 32)

    def testClone(self):
        p = self.table.props
        self.addprops(p)
        m2 = self.m.clone()
        p = m2.table(self.TABLE_NAME).props
        self.checkprops(p)

    def testAppend(self):
        self.addprops(self.table.props)
        m = msys.CreateSystem()
        m.append(self.m)
        self.checkprops(m.table(self.TABLE_NAME).props)

    def testNonbonded(self):
        t = self.m.addTable("nonbonded", 1)
        t.category = "nonbonded"
        self.m.addAtom().atomic_number = 6
        t.addTerm([self.m.atom(0)], t.params.addParam())
        self.addprops(t.props)
        with tempfile.NamedTemporaryFile(suffix=".dms") as tmp:
            msys.Save(self.m, tmp.name)
            m2 = msys.Load(tmp.name)
        self.checkprops(m2.table("nonbonded").props)

    def testDms(self):
        p = self.table.props
        mol = self.table.system
        mol.addAtom().atomic_number = 8
        mol.addAtom().atomic_number = 1
        mol.addAtom().atomic_number = 1
        param = self.table.params.addParam()
        self.table.addTerm([mol.atom(0), mol.atom(1)], param)
        self.table.addTerm([mol.atom(0), mol.atom(2)], param)

        self.addprops(p)
        self.checkprops(p)
        with tempfile.NamedTemporaryFile(suffix=".dms") as tmp:
            msys.Save(mol, tmp.name)
            mol2 = msys.Load(tmp.name)
            table = mol2.table(self.TABLE_NAME)
            self.checkprops(table.props)


class TestSpatialHash(unittest.TestCase):
    def testNearest(self):
        mol = msys.Load("tests/files/small.mae")
        pos = mol.positions.astype("f")
        box = mol.cell

        pro = mol.selectArr("fragid 0")
        wat = mol.selectArr("water")

        sph = msys.SpatialHash(pos, pro)
        new1 = sph.findNearest(1930, pos, wat)
        old1 = mol.selectArr("water and nearest 1930 to fragid 0")
        self.assertTrue((old1 == new1).all())

        sph = msys.SpatialHash(pos, pro, box=box)
        new2 = sph.findNearest(1930, pos, wat)
        old2 = mol.selectArr("water and pbnearest 1930 to fragid 0")
        self.assertTrue((old2 == new2).all())
        # make sure we really tested the pb version
        self.assertTrue((old1 != old2).any())

    def testWithin(self):
        mol = msys.Load("tests/files/small.mae")
        pos = mol.positions.astype("f")
        box = mol.cell

        pro = mol.selectArr("fragid 0")
        wat = mol.selectArr("water")

        sph = msys.SpatialHash(pos[pro])
        new1 = sph.findWithin(23.0, pos, wat)
        old1 = mol.selectArr("water and within 23.0 of fragid 0")
        self.assertTrue((old1 == new1).all())

        sph = msys.SpatialHash(pos[pro], box=box)
        new2 = sph.findWithin(23.0, pos, wat)
        old2 = mol.selectArr("water and pbwithin 23.0 of fragid 0")
        self.assertTrue((old2 == new2).all())
        self.assertFalse(old1.tolist() == old2.tolist())

    def testEmptyWithin(self):
        mol = msys.CreateSystem()
        mol.addAtom()
        self.assertEqual(mol.selectIds("within 3 of none"), [])
        self.assertEqual(mol.selectIds("pbwithin 3 of none"), [])
        self.assertEqual(mol.selectIds("exwithin 3 of none"), [])


class TestReorder(unittest.TestCase):
    def setUp(self):
        m1 = msys.CreateSystem()
        a1 = m1.addAtom()
        a2 = m1.addAtom()
        a3 = m1.addAtom()
        a4 = m1.addAtom()
        for a in m1.atoms:
            a.atomic_number = 1
        a1.addBond(a2)
        a1.addBond(a3)
        a1.addBond(a4)
        self.mol = m1

    def testIdentity(self):
        reorder.Reorder(self.mol, self.mol, "all", "all")

    def testSubset(self):
        ref = self.mol
        tgt = self.mol.clone("index 0 1 2")

        with self.assertRaises(ValueError):
            reorder.Reorder(ref, tgt, "all", "all")

        reorder.Reorder(ref, tgt, "index 0 1 2", "all")


class TestInChI(unittest.TestCase):
    gold = {
        "fused.sdf": "InChI=1/C6H8N4O/c7-6-8-3-1-2-4(11)9-5(3)10-6/h1-3,5H,(H,9,11)(H3,7,8,10)/q+1/f/h9-10H,7H2",
        "jandor.sdf": "InChI=1/C21H20N4O7S/c1-11(26)14-15-17(31-2)18(33-21-22-8-3-9-23-21)16(24(15)19(14)27)20(28)32-10-12-4-6-13(7-5-12)25(29)30/h3-9,11,14-15,17,26H,10H2,1-2H3",
    }

    gold[
        "lig.sdf"
    ] = "InChI=1/C15H27N3O2/c1-18-7-6-9-12(18)5-4-10-14-11(8-13(16-14)19-2)17-15(9,10)20-3/h9-14,16-17H,4-8H2,1-3H3/p+2/fC15H29N3O2/h16,18H/q+2"

    def testBuiltin(self):
        from msys import InChI

        for k, v in list(self.gold.items()):
            m = msys.Load("tests/files/%s" % k)
            self.assertEqual(InChI(m, SNon=True).string, v)

    def testVsites(self):
        mol = msys.Load("tests/files/tip5p.mae").clone("fragid 0")
        # for some reason our tip5p file has formal_charge=1 on everything
        for a in mol.atoms:
            a.formal_charge = 0
        self.assertEqual([a.atomic_number for a in mol.atoms], [8, 1, 1, 0, 0])
        inchi = msys.InChI(mol).string
        self.assertEqual(inchi, "InChI=1/H2O/h1H2")


class TestAmber(unittest.TestCase):
    def testStructureOnly(self):
        m1 = msys.LoadPrmTop("tests/files/molecule.prmtop")
        m2 = msys.LoadPrmTop("tests/files/molecule.prmtop", structure_only=True)
        self.assertEqual(
            [(b.first.id, b.second.id) for b in m1.bonds],
            [(b.first.id, b.second.id) for b in m2.bonds],
        )
        self.assertTrue(m1.getTable("stretch_harm") is not None)
        self.assertTrue(m2.getTable("stretch_harm") is None)

    def testWithoutTables(self):
        m = msys.Load(
            "tests/files/molecule.prmtop", structure_only=True, without_tables=False
        )
        self.assertTrue(m.getTable("stretch_harm") is not None)
        m = msys.Load(
            "tests/files/molecule.prmtop", structure_only=False, without_tables=True
        )
        self.assertTrue(m.getTable("stretch_harm") is None)

    def testHbond(self):
        m1 = msys.Load("tests/files/sys.prmtop")
        msys.ReadCrdCoordinates(m1, "tests/files/eq.rst")
        # msys.Save(m1, 'foo.dms')
        wat = m1.clone("water")
        self.assertEqual(wat.nbonds, 2 * wat.nresidues)


class TestFunkyStretch(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.mol = msys.Load("tests/files/noe.mae")

    def testFbhw(self):
        mol = self.mol
        fbhw = mol.table("stretch_fbhw")
        self.assertEqual(
            set(fbhw.params.props), set(("lower", "upper", "sigma", "beta", "fc"))
        )
        self.assertEqual(fbhw.term_props, ["group"])
        self.assertEqual([t.atoms[0].id for t in fbhw.terms], list(range(6)))
        self.assertEqual([t["group"] for t in fbhw.terms], [0] * 3 + [1] * 3)

    def testMorse(self):
        mol = self.mol
        morse = mol.table("stretch_morse")
        p = morse.params.param(0)
        self.assertEqual(p["r0"], 1.529)
        self.assertEqual(p["d"], 268)
        self.assertEqual(p["a"], 2)


class TestAlchemicalMorse(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.mol = msys.Load("tests/files/morse.mae")

    def testMorse(self):
        mol = self.mol
        morse = mol.table("alchemical_stretch_morse")
        p = morse.params.param(0)
        self.assertEqual(p["r0A"], 1.529)
        self.assertEqual(p["dA"], 268)
        self.assertEqual(p["aA"], 2)
        self.assertEqual(p["r0B"], 1.529)
        self.assertEqual(p["dB"], 0)
        self.assertEqual(p["aB"], 0)


class AtomselCoverage(unittest.TestCase):
    def setUp(self):
        """Check that all keywords are implemented"""
        m = msys.CreateSystem()
        c1 = m.addCt()
        c2 = m.addCt()
        s1 = c1.addChain()
        s2 = c2.addChain()
        r1 = s1.addResidue()
        r2 = s2.addResidue()
        a1 = r1.addAtom()
        a2 = r2.addAtom()
        a3 = r2.addAtom()
        a1.atomic_number = 6
        a2.atomic_number = 1
        a2.charge = 0.5
        a1.name = "0A"
        a2.name = "5'"
        a2.addBond(a3)
        r1.name = "A"
        r2.name = "T"

        s1.name = "A"
        s2.name = " A "
        s1.segid = "XYZ"
        s2.segid = ""
        a2.vy = 42

        m.updateFragids()
        self.mol = m

    def check(self, sel, ids):
        got = self.mol.selectIds(sel)
        self.assertEqual(got, ids, "Failed on '%s', want %s, got %s" % (sel, ids, got))

    def testEmpty(self):
        with self.assertRaises(RuntimeError):
            self.mol.select("")

    def testIons(self):
        m = msys.CreateSystem()
        a = m.addAtom()
        for resname in "Cl-", "K+":
            a.residue.name = resname
            self.assertEqual(len(m.select("legacy_ion")), 1, "failed for %s" % resname)
        a.atomic_number = 9
        self.assertEqual(m.selectIds("ion"), [0])
        for n in 0, 1, 2, 5, 6, 7, 8:
            a.atomic_number = n
        self.assertEqual(m.selectIds("ion"), [])

    def testSmallestRing(self):
        """check that r<n> conforms to spec, and test q<n> for legacy behavior"""
        m = msys.FromSmilesString("C1CC12CC(C2)")
        a = msys.AnnotatedSystem(m)

        # q<n> replaces r<n> from previous versions
        # matches "atom in ring of size n"
        p3 = msys.SmartsPattern("[C;X4;H0;q3;v4]")
        p4 = msys.SmartsPattern("[C;X4;H0;q4;v4]")
        self.assertTrue(p3.match(a))
        self.assertTrue(p4.match(a))

        # r<n> matches "atom in *smallest* ring of size n"
        p3 = msys.SmartsPattern("[C;X4;H0;r3;v4]")
        p4 = msys.SmartsPattern("[C;X4;H0;r4;v4]")
        self.assertTrue(p3.match(a))
        self.assertFalse(p4.match(a))

    def testDegree(self):
        m = msys.CreateSystem()
        a0 = m.addAtom()
        a1 = m.addAtom()
        a2 = m.addAtom()
        a0.atomic_number = 6
        a1.atomic_number = 0
        a2.atomic_number = 1
        a0.addBond(a1)
        a0.addBond(a2)
        self.assertEqual(m.selectIds("numbonds 0"), [])
        self.assertEqual(m.selectIds("numbonds 1"), [1, 2])
        self.assertEqual(m.selectIds("numbonds 2"), [0])
        self.assertEqual(m.selectIds("degree 0"), [1])
        self.assertEqual(m.selectIds("degree 1"), [0, 2])
        self.assertEqual(m.selectIds("degree 2"), [])

    def testAllBoolean(self):
        heavy = [a.id for a in self.mol.atoms if a.atomic_number > 1]
        self.check("atomicnumber > 1", heavy)
        self.check("all and atomicnumber > 1", heavy)
        self.check("atomicnumber > 1 and all", heavy)

        self.check("none and atomicnumber > 1", [])
        self.check("atomicnumber > 1 and none", [])
        self.check("none or atomicnumber > 1", heavy)
        self.check("atomicnumber > 1 or none", heavy)

    def testAnum(self):
        self.check("atomicnumber 6", [0])
        self.check("atomicnumber 1", [1])
        self.check("atomicnumber 1 6", [0, 1])
        self.check("atomicnumber 7", [])

    def testChain(self):
        self.check("chain A", [0])
        self.check("chain ' A '", [1, 2])

    def testSegid(self):
        self.check("segid XYZ", [0])
        self.check("segid ''", [1, 2])

    def testCtnumber(self):
        self.check("ctnumber 0", [])
        self.check("ctnumber 1", [0])
        self.check("ctnumber 2", [1, 2])
        self.check("ctnumber 3", [])

    def testCt(self):
        self.check("ct 0", [0])
        self.check("ct 1", [1, 2])
        self.check("ct 2", [])

    def testCharge(self):
        self.check("charge 0.0", [0, 2])
        self.check("charge 0.5", [1])
        self.check("charge==0", [0, 2])
        self.check("charge==1/2", [1])

    def testElement(self):
        self.check("element C", [0])
        self.check("element H", [1])
        self.check("element H C", [0, 1])
        self.check("element C H", [0, 1])
        self.check("element N", [])
        self.check("element XYZ", [])

    def testFragment(self):
        self.check("fragment 0", [0])
        self.check("fragment 1", [1, 2])
        self.check("fragid 0", [0])
        self.check("fragid 1", [1, 2])

    def testIndex(self):
        self.check("index 0", [0])
        self.check("index 1", [1])
        self.check("index 99", [])

    def testMass(self):
        self.check("mass 0", [0, 1, 2])
        self.check("mass<0", [])
        self.check("mass<=0", [0, 1, 2])
        self.check("mass>=0", [0, 1, 2])
        self.check("mass>0", [])

    def testName(self):
        self.check("name 0A", [0])
        self.check("name '0A'", [0])
        self.check("name 5'", [1])

    def testNumbonds(self):
        self.check("numbonds 0", [0])
        self.check("numbonds 1", [1, 2])

    def testRes(self):
        self.check("residue 0", [0])
        self.check("residue 1", [1, 2])

        self.check("resid 0", [0, 1, 2])
        self.check("resname A", [0])
        self.check("resname T", [1, 2])
        self.check('insertion ""', [0, 1, 2])

    def testPosVel(self):
        self.check("x 0", [0, 1, 2])
        self.check("y 0", [0, 1, 2])
        self.check("z 0", [0, 1, 2])
        self.check("vx 0", [0, 1, 2])
        self.check("vy 0", [0, 2])
        self.check("vz 0", [0, 1, 2])
        pos = NP.ones((3, 3), dtype="f")
        pos[2, 2] = 42
        pos[1, 1] = -3
        self.assertEqual(self.mol.selectIds("x 0", pos), [])
        self.assertEqual(self.mol.selectIds("x 1", pos), [0, 1, 2])
        self.assertEqual(self.mol.selectIds("y -3", pos), [1])
        self.assertEqual(self.mol.selectIds("z 42", pos), [2])

    def testCmp(self):
        self.check("vy == 42", [1])
        self.check("vy != 42", [0, 2])


class TestAtomsel(unittest.TestCase):
    def testCappingResiduesACE_NME(self):
        mol = msys.Load("tests/files/solvated_test.mol2")
        ace = mol.selectIds("resname ACE")
        nme = mol.selectIds("resname NME")
        tot = sorted(ace + nme)
        self.assertEqual(mol.selectIds("protein and resname ACE NME"), tot)

    def testSequence(self):
        mol = msys.Load("tests/files/2f4k.dms")
        for sel, ids in (
            ("name CA and sequence F", [73, 141, 244]),
            ("name CA and sequence RS", [199, 223]),
            ("name CA and sequence R..F", [199, 223, 234, 244]),
            ("name CA and sequence FG..R.A", [141, 161, 168, 185, 199, 223, 234]),
            ("name CA and sequence F K", [73, 93, 141, 244, 473, 510]),
        ):
            self.assertEqual(mol.selectIds(sel), ids)

    def testSmarts(self):
        mol = msys.Load("tests/files/ch4.dms")
        self.assertEqual(mol.selectIds("smarts C"), [0])
        self.assertEqual(mol.selectIds("smarts '[H]'"), [1, 2, 3, 4])

    def testPrecedence(self):
        mol = msys.Load("tests/files/2f4k.dms")
        sels = {
            "within 3 of name CA or name CB": 482,
            "within 3 of resname ALA and noh": 63,
            "not resname ALA or hydrogen": 13198,
            "not resname ALA and hydrogen": 8697,
            "nearest 5 to name CA or name CB": 5,
        }
        for sel, cnt in sels.items():
            self.assertEqual(len(mol.selectIds(sel)), cnt)

    # def testParamtype(self):
    # mol=msys.Load('tests/files/2f4k.dms')
    # self.assertEqual(mol.selectIds('paramtype nonbonded'), [])

    def testTrajectory(self):
        mol = msys.Load("tests/files/alanin.pdb")
        trj = molfile.dcd.read("tests/files/alanin.dcd")
        for f in trj.frames():
            box = NP.diag(NP.max(f.pos, 0) - NP.min(f.pos, 0))
            id1 = mol.selectIds("within 5 of residue 3", pos=f.pos, box=box)
            id2 = mol.selectIds("pbwithin 7 of residue 3", pos=f.pos, box=box)

            mol.setPositions(f.pos)
            mol.setCell(box)
            id3 = mol.selectIds("within 5 of residue 3")
            id4 = mol.selectIds("pbwithin 7 of residue 3")

            self.assertEqual(id1, id3)
            self.assertEqual(id2, id4)

    def testBaselines(self):
        with open("tests/atomsel_tests.json") as f:
            d = json.load(f)
        for p, base in list(d.items()):
            mol = msys.Load(str(p))
            for sel, old in base:
                sel = str(sel)
                new = mol.selectIds(sel)
                self.assertEqual(
                    old,
                    new,
                    "failed on '%s': oldlen %d newlen %d" % (sel, len(old), len(new)),
                )


class TestSdf(unittest.TestCase):
    def testFailedOpenWrite(self):
        mol = msys.Load("tests/files/isotope.sdf")
        with self.assertRaises(RuntimeError):
            mol.save("/this/will/not/work.sdf")

    def testIsotopes(self):
        mol = msys.Load("tests/files/isotope.sdf")
        assert "isotope" in mol.atom_props
        for i in 1, 2:
            assert mol.atom(i)["isotope"] == 2

    def testZeroBondOrder(self):
        mol = msys.CreateSystem()
        c = mol.addAtom()
        c.atomic_number = 17
        c.formal_charge = -1
        v = mol.addAtom()
        v.atomic_number = 0
        bnd = c.addBond(v)
        bnd.order = 0
        sdf = msys.FormatSDF(mol)
        new = next(msys.ParseSDF(sdf))
        self.assertEqual(new.bond(0).order, 0)
        self.assertEqual(new.atom(0).formal_charge, -1)

    def testIter(self):
        with open("tests/files/cofactors.sdf") as fp:
            s = fp.read()
        for m in msys.ParseSDF(s):
            pass

    def testCoordPrecision(self):
        tmp = tempfile.NamedTemporaryFile(suffix=".sdf")
        mol = msys.CreateSystem()
        for i in range(999):
            a = mol.addAtom()
            a.atomic_number = 6
            a.x = 0.0001 * (i + 1)
        msys.Save(mol, tmp.name)
        s = msys.FormatSDF(mol)
        lines = s.split("\n")[4 : 999 + 4]
        xs = [line.split()[0] for line in lines]
        mol2 = msys.Load(tmp.name)
        for i in range(999):
            self.assertEqual(xs[i], "%.4f" % ((i + 1) / 10000.0))
            a = mol.atom(i)
            self.assertEqual(a.x, 0.0001 * (i + 1))

    def testExtraLineSdfLoad(self):
        path = "tests/files/extra_line_m_end.sdf"
        for x in msys.LoadMany(path):
            self.assertFalse(x is None)

    def testExtraLineSdfScan(self):
        path = "tests/files/extra_line_m_end.sdf"
        for x in msys.LoadMany(path):
            self.assertFalse(x is None)

    def testFormalCharge(self):
        mol = msys.Load("tests/files/lig.sdf")
        self.assertEqual(mol.atom(12).formal_charge, 1)
        self.assertEqual(mol.atom(19).formal_charge, 1)
        self.assertEqual(mol.atom(20).formal_charge, 0)
        self.assertEqual(mol.atom(21).formal_charge, 0)

    def testSingleCt(self):
        tmp = tempfile.NamedTemporaryFile(suffix=".sdf")
        path = tmp.name
        mol = msys.CreateSystem()
        a = mol.addAtom()
        mol.ct(0).name = "XYZ"
        a.atomic_number = 11
        msys.Save(mol, path)
        mol2 = msys.Load(path)
        self.assertEqual(mol2.ct(0).name, "XYZ")

    def testFunnyDoubleLikeThing(self):
        mol = msys.Load("tests/files/multiline.sdf")
        self.assertEqual(mol.ct(0)["PUBCHEM_CONFORMER_ID"], "0001E10300000001")

    def testMultiline(self):
        mol = msys.Load("tests/files/multiline.sdf")
        self.assertEqual(
            mol.ct(0)["PUBCHEM_SHAPE_FINGERPRINT"],
            "21015797 1 18410576218592473837\n260 1 18410856563971680135",
        )

    def testMultilineWithBlankOrSpace(self):
        mol = msys.Load("tests/files/multiline.sdf")
        self.assertEqual(
            mol.ct(0)["PUBCHEM_SHAPE_MULTIPOLES"],
            "\n".join(
                map(
                    str,
                    [36.18, 0.98, 0.57, 0.57, 0.03, "", 0, " ", 0, 0.03, 0, 0.03, 0, 0],
                )
            ),
        )

    def testMultipleCt(self):
        tmp = tempfile.NamedTemporaryFile(suffix=".sdf")
        path = tmp.name
        ct = msys.CreateSystem()
        a = ct.addAtom()
        a.atomic_number = 11
        mol = msys.CreateSystem()
        mol.append(ct)
        mol.append(ct)
        mol.ct(0).name = "XYZ"
        mol.ct(1).name = "ABC"
        self.assertEqual(mol.ct(0).name, "XYZ")
        self.assertEqual(mol.ct(1).name, "ABC")

        msys.Save(mol, path)
        mol2 = msys.Load(path)
        self.assertEqual(mol2.ct(0).name, "XYZ")
        self.assertEqual(mol2.ct(1).name, "ABC")
        for i, mol in enumerate(msys.LoadMany(path)):
            self.assertEqual(mol.ct(0).name, ("XYZ", "ABC")[i])
            self.assertEqual(mol.name, ("XYZ", "ABC")[i])

    def testNoDelim(self):
        msys.Load("tests/files/no-delim.sdf")

    def testNoDelim2(self):
        mol = msys.Load("tests/files/no-delim2.sdf")
        self.assertEqual(mol.ct(0)["Name"], "dUMP anion\nanother line")


class TestMolfile(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        r = molfile.pdb.read("tests/files/3RYZ.pdb")
        cls.atoms = r.atoms
        cls.frame = r.frame(0)

    def testOccupancy(self):
        ids = [i for i, a in enumerate(self.atoms) if a.occupancy < 0.3]
        self.assertEqual(ids, [490, 494, 496, 498, 500, 502, 504])

    def testBfactor(self):
        ref = [2182, 2191, 2208, 2211, 2223, 2259, 2264, 2279, 2393, 2400, 2441]
        ids = [i for i, a in enumerate(self.atoms) if a.bfactor > 50]
        self.assertEqual(ids, ref)
        with tempfile.NamedTemporaryFile(suffix=".pdb") as tmp:
            molfile.pdb.write(tmp.name, self.atoms).frame(self.frame).close()
            r = molfile.pdb.read(tmp.name)
            atoms = r.atoms
            ids = [i for i, a in enumerate(atoms) if a.bfactor > 50]
            self.assertEqual(ids, ref)

    def testAltloc(self):
        ids = [i for i, a in enumerate(self.atoms) if a.altloc == "A"]
        self.assertEqual(
            ids,
            [
                161,
                165,
                167,
                169,
                171,
                489,
                493,
                495,
                497,
                499,
                501,
                503,
                1448,
                1452,
                1454,
                1456,
                1458,
                1460,
                1698,
                1702,
                1704,
            ],
        )

    def testCell(self):
        self.assertEqual(
            list(self.frame.box.diagonal()), [42.408, 41.697, 69.9599803728577]
        )


class TestPdb(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.mol = msys.Load("tests/files/3RYZ.pdb")

    def testOccupancy(self):
        ids = self.mol.selectIds("occupancy < 0.3")
        self.assertEqual(ids, [490, 494, 496, 498, 500, 502, 504])

    def testBfactor(self):
        ref = [2182, 2191, 2208, 2211, 2223, 2259, 2264, 2279, 2393, 2400, 2441]
        ids = self.mol.selectIds("bfactor > 50")
        self.assertEqual(ids, ref)
        with tempfile.NamedTemporaryFile(suffix=".pdb") as tmp:
            msys.Save(self.mol, tmp.name)
            mol = msys.Load(tmp.name)
            ids = mol.selectIds("bfactor > 50")
            self.assertEqual(ids, ref)

    def testAltloc(self):
        self.assertEqual(
            self.mol.selectIds("altloc A"),
            [
                161,
                165,
                167,
                169,
                171,
                489,
                493,
                495,
                497,
                499,
                501,
                503,
                1448,
                1452,
                1454,
                1456,
                1458,
                1460,
                1698,
                1702,
                1704,
            ],
        )

    def testSave(self):
        with tmpfile(suffix=".pdb") as tmp:
            msys.Save(self.mol, tmp.name)
            oldmol = self.mol
            self.mol = msys.Load(tmp.name)
            try:
                self.testAltloc()
                self.testBfactor()
                self.testOccupancy()
            finally:
                self.mol = oldmol

    def testTer(self):
        mol = msys.CreateSystem()
        c1 = mol.addChain()
        c2 = mol.addChain()
        for c in mol.chains:
            c.name = "A"
            r = c.addResidue()
            r.name = "ALA"
            a = r.addAtom()
            a.name = "CA"
            a.atomic_number = 6
        tmp = tempfile.NamedTemporaryFile(suffix=".pdb")
        path = tmp.name
        msys.SavePDB(mol, path, reorder=True)
        m2 = msys.Load(path)
        self.assertEqual(m2.nchains, 2)

    def testReorder(self):
        mol = msys.Load("tests/files/order.mol2")
        tmp = tempfile.NamedTemporaryFile(suffix=".pdb")
        path = tmp.name

        msys.SavePDB(mol, path, reorder=True)
        pdb = msys.Load(path)
        assert [a.name for a in mol.atoms] != [a.name for a in pdb.atoms]

        msys.SavePDB(mol, path, reorder=False)
        pdb = msys.Load(path)
        assert [a.name for a in mol.atoms] == [a.name for a in pdb.atoms]

        msys.Save(mol, path)
        pdb = msys.Load(path)
        assert [a.name for a in mol.atoms] == [a.name for a in pdb.atoms]

    def testCell(self):
        self.assertEqual(
            list(self.mol.cell.diagonal()), [42.408, 41.697, 69.9599803728577]
        )


class TestValidate(unittest.TestCase):
    def testNoKnot(self):
        mol = msys.Load("tests/files/jandor.sdf")
        knot.FindKnots(mol, verbose=False)

    def testPeriodicKnot(self):
        mol = msys.Load("tests/files/knot.mae")
        results = knot.FindKnots(mol, verbose=False)
        self.assertEqual(len(results), 2)

    def testExcludedKnot(self):
        mol = msys.Load("tests/files/excluded_knot.dms")
        without_ignoring = knot.FindKnots(mol, verbose=False)
        self.assertEqual(len(without_ignoring), 1)
        with_ignoring = knot.FindKnots(mol, ignore_excluded_knots=True, verbose=False)
        self.assertEqual(len(with_ignoring), 0)

    def testUntieKnot(self):
        mol = msys.Load("tests/files/knot.mae")
        results = knot.FindKnots(mol, verbose=False)
        self.assertEqual(len(results), 2)
        success = knot.UntieKnots(mol, verbose=False)
        results_after_untying = knot.FindKnots(mol, verbose=False)
        self.assertEqual(len(results_after_untying), 0)
        self.assertTrue(success)


class Tools(unittest.TestCase):
    @staticmethod
    def make_systems():
        tmp1 = tmpfile(suffix=".dms")
        tmp2 = tmpfile(suffix=".dms")
        mol1 = msys.CreateSystem()
        mol2 = msys.CreateSystem()
        mol1.addAtom().atomic_number = 1
        mol2.addAtom().atomic_number = 2
        msys.Save(mol1, tmp1.name)
        msys.Save(mol2, tmp2.name)
        return tmp1, tmp2

    def testCell(self):
        self.assertEqual(
            list(self.mol.cell.diagonal()), [42.408, 41.697, 69.9599803728577]
        )


class TestValidate(unittest.TestCase):
    def testNoKnot(self):
        mol = msys.Load("tests/files/jandor.sdf")
        knot.FindKnots(mol, verbose=False)

    def testPeriodicKnot(self):
        mol = msys.Load("tests/files/knot.mae")
        results = knot.FindKnots(mol, verbose=False)
        self.assertEqual(len(results), 2)

    def testExcludedKnot(self):
        mol = msys.Load("tests/files/excluded_knot.dms")
        without_ignoring = knot.FindKnots(mol, verbose=False)
        self.assertEqual(len(without_ignoring), 1)
        with_ignoring = knot.FindKnots(mol, ignore_excluded_knots=True, verbose=False)
        self.assertEqual(len(with_ignoring), 0)

    def testUntieKnot(self):
        mol = msys.Load("tests/files/knot.mae")
        results = knot.FindKnots(mol, verbose=False)
        self.assertEqual(len(results), 2)
        success = knot.UntieKnots(mol, verbose=False)
        results_after_untying = knot.FindKnots(mol, verbose=False)
        self.assertEqual(len(results_after_untying), 0)
        self.assertTrue(success)


class Tools(unittest.TestCase):
    @staticmethod
    def make_systems():
        tmp1 = tmpfile(suffix=".dms")
        tmp2 = tmpfile(suffix=".dms")
        mol1 = msys.CreateSystem()
        mol2 = msys.CreateSystem()
        mol1.addAtom().atomic_number = 1
        mol2.addAtom().atomic_number = 2
        msys.Save(mol1, tmp1.name)
        msys.Save(mol2, tmp2.name)
        return tmp1, tmp2

    def testDiff(self):
        from msys import dump

        file1, file2 = self.make_systems()
        self.assertEqual(dump.diff_main([file1.name, file1.name]), 0)
        self.assertEqual(dump.diff_main([file2.name, file1.name], quiet=True), 1)
        self.assertEqual(dump.diff_main([file1.name, file2.name], quiet=True), 1)

    def testDmsDiff(self):
        file1, file2 = self.make_systems()
        self.assertEqual(
            subprocess.call(
                ["dms-diff", file1.name, file1.name], stdout=subprocess.PIPE
            ),
            0,
        )
        self.assertEqual(
            subprocess.call(
                ["dms-diff", file2.name, file1.name], stdout=subprocess.PIPE
            ),
            1,
        )
        self.assertEqual(
            subprocess.call(
                ["dms-diff", file1.name, file2.name], stdout=subprocess.PIPE
            ),
            1,
        )

    def testDmsSelect(self):
        file1, file2 = self.make_systems()
        self.assertEqual(
            subprocess.call(
                ["dms-select", file1.name, "-o", file2.name], stderr=subprocess.PIPE
            ),
            0,
        )


class Main(unittest.TestCase):
    def testInfo(self):
        from msys.info import print_info
        mol = msys.Load("tests/files/ww.dms")
        print_info(mol)


    def testMixedCtPropertyTypes(self):
        mol = msys.CreateSystem()
        ct0 = mol.addCt()
        ct1 = mol.addCt()
        ct2 = mol.addCt()
        ct3 = mol.addCt()
        ct0["x"] = 1
        ct1["x"] = 3.5
        ct2["x"] = "1"
        ct3["x"] = "x"

        s = msys.FormatDMS(mol)
        m = msys.LoadDMS(buffer=s)

        assert m.ct(0)["x"] == 1
        assert m.ct(1)["x"] == 3.5
        assert m.ct(2)["x"] == 1  # sadly, numeric conversion occurs
        assert m.ct(3)["x"] == "x"

    def testSave(self):
        mol1 = msys.Load("tests/files/ch4.dms")
        tmp = tmpfile(suffix=".dms")
        mol1.save(mol1, tmp.name)
        mol2 = msys.Load(tmp.name)
        assert mol1.hash() == mol2.hash()

    def testMol2NoResidues(self):
        mol = msys.Load("tests/files/m_112.mol2")
        self.assertEqual(mol.nresidues, 1)
        self.assertEqual(mol.natoms, 18)

    def testMol2NonconsecutiveResid(self):
        mol = msys.Load("tests/files/nonconsecutive-resid.mol2")
        self.assertEqual(mol.nresidues, 301)

    def testMol2ChainSplit(self):
        mol = msys.Load("tests/files/complex_aligned.mol2")
        assert mol.nchains == 2, "got %d chains, expected 2" % mol.nchains

    def testMol2ChainsMerge(self):
        mol = msys.Load("tests/files/solvated_test.mol2")
        assert mol.nchains == 6, "got %d chains, expected 6" % mol.nchains
        assert mol.ncts == 2
        with tmpfile(suffix=".mol2") as tmp:
            msys.Save(mol, tmp.name)
            mol2 = msys.Load(tmp.name)
        assert mol2.nchains == 6, "got %d chains, expected 6" % mol2.nchains
        assert mol2.ncts == 2, "got %d cts, expected 2" % mol2.ncts

    def testAddProvenance(self):
        mol = msys.CreateSystem()
        prov = msys._msys.Provenance()
        attrs = "version", "timestamp", "user", "workdir", "cmdline", "executable"
        for attr in attrs:
            setattr(prov, attr, "the %s" % attr)
        mol._ptr.addProvenance(prov)
        with tempfile.NamedTemporaryFile(suffix=".dms") as tmp:
            msys.Save(mol, tmp.name)
            new = msys.Load(tmp.name)
        self.assertEqual(len(new.provenance), 2)
        for attr in attrs:
            self.assertEqual(getattr(prov, attr), getattr(new.provenance[0], attr))

    def testReplaceProvenance(self):
        mol = msys.Load("tests/files/pseudo.dms")
        prov = mol.provenance
        prov = prov[1:3]
        assert len(prov) == 2
        mol.provenance = prov
        attrs = "version", "timestamp", "user", "workdir", "cmdline", "executable"
        assert len(prov) == len(mol.provenance)
        for old, new in zip(prov, mol.provenance):
            for attr in attrs:
                assert getattr(old, attr) == getattr(new, attr)

    def testCloneWholeFragments(self):
        mol = msys.FromSmilesString("CC")
        mol.clone()
        mol.clone(forbid_broken_bonds=True)
        for i in range(mol.natoms):
            with self.assertRaises(msys.BrokenBondsError):
                mol.clone([i], forbid_broken_bonds=True)
            mol.clone([i])
            mol.clone([i], forbid_broken_bonds=False)

    def testFindDistinctFragments(self):
        import random

        random.seed(99)
        frag1 = msys.FromSmilesString("CCCC")
        frag2 = msys.FromSmilesString("CC(C)C")
        atoms = frag1.atoms
        random.shuffle(atoms)
        frag1b = frag1.clone(atoms)

        assert msys.InChI(frag1).string != msys.InChI(frag2).string
        assert msys.Graph(frag1).hash() == msys.Graph(frag2).hash()
        assert [a.atomic_number for a in frag1.atoms] != [
            a.atomic_number for a in frag1b.atoms
        ]

        mol = msys.CreateSystem()
        mol.append(frag2)
        mol.append(frag1)
        mol.append(frag1b)
        mol.append(frag2)
        mol.append(frag1)

        result = msys.FindDistinctFragments(mol)
        self.assertEqual(result, {0: [0, 3], 1: [1, 2, 4]})

    def testFindDistinctFragmentsStereo(self):
        mol = msys.Load("tests/files/stereo.sdf")
        frags = msys.FindDistinctFragments(mol)
        self.assertEqual(len(frags), 1)

        frags = msys.FindDistinctFragments(mol, key="inchi")
        self.assertEqual(len(frags), 2)

    def testSortAtoms(self):
        mol = msys.Load("tests/files/2f4k.dms")
        atoms = mol.atoms
        random.shuffle(atoms)
        ids = [a.id for a in atoms]
        ref = sorted(ids)
        self.assertNotEqual(ids, ref)
        atoms.sort()
        ids = [a.id for a in atoms]
        self.assertEqual(ids, ref)

    def testSortBonds(self):
        mol = msys.Load("tests/files/2f4k.dms")
        bonds = mol.bonds
        random.shuffle(bonds)
        ids = [(b.first.id, b.second.id) for b in bonds]
        ref = sorted(ids)
        self.assertNotEqual(ids, ref)
        bonds.sort()
        ids = [(b.first.id, b.second.id) for b in bonds]
        self.assertEqual(ids, ref)

    def testSortMismatch(self):
        import itertools

        mol = msys.Load("tests/files/ch4.dms")
        atm = mol.atom(0)
        bnd = mol.bond(0)
        for i, j in itertools.permutations((mol, atm, bnd), 2):
            with self.assertRaises(TypeError):
                [i, j].sort()
                [j, i].sort()

    def testCapsule(self):
        mol = msys.Load("tests/files/2f4k.dms")
        ptr = mol._ptr
        cap = msys._msys.SystemPtr.asCapsule(ptr)
        ptr2 = msys._msys.SystemPtr.fromCapsule(cap)
        self.assertTrue(ptr is ptr2)
        self.assertEqual(ptr, ptr2)

    def testCapsule2(self):
        mol = msys.Load("tests/files/2f4k.dms")
        cap = mol.asCapsule()
        mol2 = msys.System.fromCapsule(cap)
        self.assertFalse(mol is mol2)
        self.assertEqual(mol, mol2)

    def testPickle(self):
        import pickle as pkl

        old = msys.Load("tests/files/2f4k.dms")
        s = pkl.dumps(old, pkl.HIGHEST_PROTOCOL)
        new = pkl.loads(s)
        self.assertEqual(old.natoms, new.natoms)
        self.assertEqual(old.nbonds, new.nbonds)
        self.assertEqual(old.table_names, new.table_names)

        # can pickle residues
        oldres = old.residues[:10]
        s = pkl.dumps(oldres, pkl.HIGHEST_PROTOCOL)
        newres = pkl.loads(s)
        self.assertEqual([a.resid for a in oldres], [a.resid for a in newres])

    def testPickleProvenance(self):
        import copy

        old = msys.CreateSystem()
        new = copy.deepcopy(old)
        self.assertEqual(new.provenance, [])

    def testFormatDMS(self):
        old = msys.Load("tests/files/2f4k.dms")
        dms = msys.FormatDMS(old)
        new = msys.LoadDMS(buffer=dms)
        self.assertEqual(old.hash(), new.hash())

    def testGuessBonds(self):
        mol = msys.Load("tests/files/2f4k.dms")
        mol.guessBonds()

    def testGeometry(self):
        p = NP.array(((1, 0, 0), (2, 3, 1), (3, 1, 0), (2, 2, 2)), "d")
        a, b, c, d = p
        rad = msys.CalcDistance(c, d)
        theta = msys.CalcAngle(b, c, d)
        phi = msys.CalcDihedral(a, b, c, d)
        tst = msys.ApplyDihedralGeometry(a, b, c, rad, theta, phi)
        NP.testing.assert_almost_equal(tst, d)

    def testMaePrecision(self):
        """ensure mae doesn't truncate coordinate precision"""
        mol = msys.CreateSystem()
        atm = mol.addAtom()
        atm.atomic_number = 17
        pos = [1.0000001, 20.447919, -550.12999]
        atm.pos = pos
        assert atm.pos.tolist() == pos
        dms = msys.LoadDMS(buffer=msys.FormatDMS(mol))
        assert dms.atom(0).pos.tolist() == pos
        mae = msys.LoadMAE(buffer=msys.SerializeMAE(mol))
        assert mae.atom(0).pos.tolist() == pos


    def testMaeNoncontiguous(self):
        """disallow writing mae when it would change atom order"""
        m = msys.CreateSystem()
        ct0 = m.addCt()
        ct1 = m.addCt()
        a0 = m.addAtom()
        a1 = m.addAtom()
        a2 = m.addAtom()
        a0.name = "a0"
        a1.name = "a1"
        a2.name = "a2"

        self.assertEqual(ct0.natoms, 3)
        self.assertEqual(ct1.natoms, 0)
        msys.SerializeMAE(m)

        a1.residue.chain.ct = ct1
        self.assertEqual(ct0.natoms, 2)
        self.assertEqual(ct1.natoms, 1)
        with self.assertRaises(RuntimeError):
            msys.SerializeMAE(m)

        msys.SerializeMAE(m, allow_reorder_atoms=True)

    def testMaeBonds(self):
        m = msys.CreateSystem()
        o = m.addAtom()
        h = m.addAtom()
        p = m.addAtom()
        o.addBond(h)
        o.addBond(p)
        o.atomic_number = 8
        h.atomic_number = 0
        p.atomic_number = 1
        tmp = tempfile.NamedTemporaryFile(suffix=".mae")
        path = tmp.name
        msys.Save(m, path)
        mol = msys.Load(path)
        self.assertEqual(mol.nbonds, 1)
        self.assertEqual(mol.bond(0).first.id, 0)
        self.assertEqual(mol.bond(0).second.id, 1)
        self.assertEqual(mol.atom(2).atomic_number, 0)

    def testTableProps(self):
        m = msys.CreateSystem()
        a = m.addAtom()
        t = m.addTable("foo", 1)
        t.category = "nonbonded"
        p = t.params.addParam()
        t.addTerm([a], p)
        t.props["foo"] = 1.0
        t.props["bar"] = 1.0

        with tempfile.NamedTemporaryFile(suffix=".dms") as fp:
            dst = fp.name
            msys.Save(m, dst)
            m2 = msys.Load(dst)
            msys.Save(m2, dst)

    def testSave(self):
        import operator as op

        old = msys.Load("tests/files/jandor.sdf")
        aget = op.attrgetter("name", "atomic_number", "charge")
        bget = op.attrgetter("first.id", "second.id", "order")
        for fmt in "dms mae mol2 sdf".split():
            tmp = tempfile.NamedTemporaryFile(suffix=".%s" % fmt)
            dst = tmp.name
            msys.Save(old, dst)
            new = msys.Load(dst)
            self.assertEqual(list(map(aget, old.atoms)), list(map(aget, new.atoms)))
            self.assertEqual(list(map(bget, old.bonds)), list(map(bget, new.bonds)))

    def testSaveAppend(self):
        import operator as op

        old = msys.Load("tests/files/jandor.sdf")
        aget = op.attrgetter("name", "atomic_number", "charge")
        bget = op.attrgetter("first.id", "second.id", "order")
        for fmt in "dms mae mol2 sdf".split():
            tmp = tempfile.NamedTemporaryFile(suffix=".%s" % fmt)
            dst = tmp.name
            msys.Save(old, dst)
            msys.Save(old, dst, append=True)
            new = msys.Load(dst)
            self.assertEqual(new.ncts, 2)
            for ct in new.cts:
                ct = new.clone(a.id for a in ct.atoms)
                self.assertEqual(list(map(aget, old.atoms)), list(map(aget, ct.atoms)))
                self.assertEqual(list(map(bget, old.bonds)), list(map(bget, ct.bonds)))

    def testBadId(self):
        m = msys.CreateSystem()
        with self.assertRaises(ValueError):
            m.chain("A")
        with self.assertRaises(ValueError):
            m.residue("A")
        with self.assertRaises(TypeError):
            m.atom("A")

    def testTermSystem(self):
        m = msys.CreateSystem()
        m.addAtom()
        s = m.addTableFromSchema("stretch_harm")
        s.addTerm((m.atom(0), m.atom(0)))
        self.assertEqual(m, s.term(0).system)

    def testRigidExplicit(self):
        m = msys.CreateSystem()
        nmax = 9
        for i in range(nmax):
            m.addAtom()
        for i in range(2, nmax + 1):
            s = m.addTableFromSchema("rigid_explicit%d" % i)
            self.assertEqual(s.natoms, i)
            self.assertEqual(s.params.nprops, 3 * i)

    def testLoadByJobid(self):
        x = msys.Load(6181535)
        self.assertEqual(x.natoms, 12960)

    def testParamSystem(self):
        m = msys.CreateSystem()
        m.addAtom()
        s = m.addTableFromSchema("stretch_harm")
        p = s.params.addParam()
        with self.assertRaises(AttributeError):
            p.system

    def testCtOrder(self):
        """ atom order in ct should be the same as the original order """
        mol = msys.CreateSystem()
        r1 = mol.addResidue()
        r2 = mol.addResidue()
        a1 = r1.addAtom()
        a2 = r2.addAtom()
        a3 = r1.addAtom()
        a4 = r2.addAtom()
        self.assertEqual(mol.ct(0).atoms, mol.atoms)

    def testCt(self):
        m = msys.CreateSystem()
        self.assertEqual(m.cts, [])
        self.assertEqual(m.ncts, 0)

        ct1 = m.addCt()
        self.assertEqual(ct1.id, 0)
        self.assertEqual(m.cts, [ct1])
        self.assertEqual(m.ncts, 1)

        ct2 = m.addCt()
        self.assertEqual(ct2.id, 1)
        self.assertEqual(m.cts, [ct1, ct2])
        self.assertEqual(m.ncts, 2)

        ct1.remove()
        self.assertEqual(m.cts, [ct2])
        self.assertEqual(m.ncts, 1)

        self.assertEqual(ct2.name, "")
        ct2.name = "jrg"
        self.assertEqual(ct2.name, "jrg")

        # newly added chains get added to the first ct
        chn = m.addChain()
        self.assertEqual(chn.ct, ct2)

        # newly added atoms get added to the first ct
        atm = m.addAtom()
        self.assertEqual(atm.residue.chain.ct, ct2)

        ct2.remove()
        chn = m.addChain()
        ct3 = m.addCt()
        self.assertEqual(chn.ct, m.cts[0])
        self.assertEqual(m.cts, [chn.ct, ct3])

        ct3["foobar"] = 32
        self.assertEqual(ct3["foobar"], 32)
        del ct3["foobar"]
        with self.assertRaises(KeyError):
            ct3["foobar"]

    def testCtAppend(self):
        m = msys.CreateSystem()
        m.addCt().name = "a"
        m2 = msys.CreateSystem()
        ct = m2.addCt()
        ct.name = "b"
        ct.addChain().name = "c"
        m3 = msys.CreateSystem()
        m4 = msys.CreateSystem()
        m4.addCt().name = "d"

        m.append(m2)
        m.append(m3)
        m.append(m4)
        self.assertEqual(m.ncts, 3)
        self.assertEqual([c.name for c in m.cts], ["a", "b", "d"])
        self.assertEqual(m.ct(1).chains[0].name, "c")

    def testAppendToCt(self):
        m = msys.CreateSystem()
        m.addCt().name = "a"
        m.addCt().name = "d"
        m.ct(1).addChain().name = "x"

        m2 = msys.CreateSystem()
        m2.cell[0][0] = 42
        ct = m2.addCt()
        ct.name = "b"
        ct.addChain().name = "c"

        m.ct(1).append(m2)
        self.assertEqual(m.cell[0][0], 0)
        self.assertEqual(m.ncts, 2)
        self.assertEqual(m.ct(0).name, "a")
        self.assertEqual(m.ct(1).name, "d")
        self.assertEqual(m.ct(1).nchains, 2)
        self.assertEqual(m.ct(1).chains[0].name, "x")
        self.assertEqual(m.ct(1).chains[1].name, "c")

    def testCtClone(self):
        m1 = msys.CreateSystem()
        m2 = msys.CreateSystem()
        m3 = msys.CreateSystem()
        for i in range(4):
            m1.addAtom().name = "%d" % i
            m2.addAtom().name = "%d" % (i + 4)
            m3.addAtom().name = "%d" % (i + 8)
        m1.append(m2)
        m1.append(m3)
        m = m1.clone()
        self.assertEqual(m1.ncts, 3)
        self.assertEqual(m.ncts, 3)
        m = m1.clone("index 0 1 9")
        self.assertEqual(m.ncts, 2)
        self.assertEqual([x.name for x in m.atoms], ["0", "1", "9"])

        # merge cts
        self.assertEqual(m.ct(0).nchains, 2)
        self.assertEqual(m.ct(1).nchains, 1)
        for chn in m.ct(0).chains:
            chn.ct = m.ct(1)
        self.assertEqual(m.ct(0).nchains, 0)
        self.assertEqual(m.ct(1).nchains, 3)
        m = m.clone()
        self.assertEqual(m.ncts, 1)

    def testMaeM2ioVersionAtEnd(self):
        msys.LoadMany("tests/files/m2io_at_end.mae")

    def testMaeMany(self):
        ct1 = msys.CreateSystem()
        ct2 = msys.CreateSystem()
        ct1.addAtom().atomic_number = 1
        ct2.addAtom().atomic_number = 2
        ct2.addAtom().atomic_number = 3
        tmp = tempfile.NamedTemporaryFile(suffix=".mae")
        dst = tmp.name
        msys.SaveMAE([ct1, ct2], dst)
        cts = [x for x in msys.LoadMany(dst)]
        self.assertEqual(len(cts), 2)
        self.assertEqual(cts[0].natoms, 1)
        self.assertEqual(cts[1].natoms, 2)

    def testPdbMany(self):
        l_sys = list(msys.LoadMany("tests/files/1DUF.pdb"))
        self.assertTrue(len(l_sys) == 5)
        sz = l_sys[0].positions.shape[0]
        self.assertTrue(sz > 0)
        for s in l_sys:
            self.assertTrue(s.positions.shape[0] == sz)

    def testMultilineCtProp(self):
        m = msys.CreateSystem()
        m.addAtom().atomic_number = 1
        prop = """
        abc
        123
        """
        m.ct(0)["prop"] = prop
        tmp = tempfile.NamedTemporaryFile(suffix=".mae")
        dst = tmp.name
        msys.SaveMAE(m, dst)
        m = msys.Load(dst)
        self.assertEqual(m.ct(0)["prop"], prop)

    def testBadAtomsel(self):
        mol = msys.CreateSystem()
        mol.addAtom()
        for s in "residue X", "residue 3.0", "residue 3.5", "mass a":
            with self.assertRaises(RuntimeError):
                mol.select(s)
        with self.assertRaises(RuntimeError):
            mol.select("residue 999999999999999999999999999")
        # with self.assertRaises(RuntimeError):
        # mol.select('residue -99999999')

    def testAtomsel(self):
        ww = "tests/files/ww.dms"
        mol = msys.Load(ww)
        ref = mol.atomsel("residue 1")
        sel = mol.atomsel("residue 3")
        oldrms = 6.0121471906466164
        newrms = 0.85438972060921559
        pos = sel.getPositions()
        self.assertAlmostEqual(ref.currentRMSD(sel), oldrms)
        self.assertAlmostEqual(ref.currentRMSD(pos), oldrms)
        self.assertAlmostEqual(ref.alignedRMSD(sel), newrms)

        ref.alignCoordinates(pos)
        self.assertAlmostEqual(ref.currentRMSD(pos), newrms)
        self.assertAlmostEqual(ref.currentRMSD(sel), oldrms)
        self.assertAlmostEqual(ref.alignedRMSD(pos), newrms)

        mol2 = msys.Load(ww)
        sel2 = mol2.atomsel("residue 3")
        self.assertAlmostEqual(ref.currentRMSD(sel2), oldrms)
        ref.alignCoordinates(sel2)
        self.assertAlmostEqual(ref.currentRMSD(sel2), newrms)

    def testAtomselAsList(self):
        ww = "tests/files/ww.dms"
        mol = msys.Load(ww)
        self.assertEqual(mol.atomsel("index 27 25").ids.tolist(), [25, 27])
        self.assertEqual(mol.atomsel((25, 27)).ids.tolist(), [25, 27])
        self.assertEqual(mol.atomsel([]).ids.tolist(), [])
        self.assertEqual(mol.atomsel(None).ids.tolist(), [])

    def testEmptyAtomsel(self):
        A = msys.Load("tests/files/alanin.pdb")
        B = A.clone()
        with self.assertRaises(RuntimeError):
            A.atomsel("none").alignCoordinates(B.atomsel("none"))

    def testRadius(self):
        for i, r in (1, 1.1), (6, 1.7), (19, 1.76):
            self.assertEqual(msys.RadiusForElement(i), r)

    def testMass(self):
        for i, m in (-1, 0), (0, 0), (1, 1.00794), (19, 39.0983):
            self.assertEqual(msys.MassForElement(i), m)
            self.assertEqual(msys.GuessAtomicNumber(m), max(i, 0))

    def testElement(self):
        for i, a, n, e in (
            (1, "H", "H", 2.3),
            (2, "He", "He", 4.16),
            (19, "K", "K", 0.734),
            (6, "C1", "C", 2.544),
            (17, "Cl2", "Cl", 2.869),
            (20, "Ca3", "Ca", 1.034),
        ):
            self.assertEqual(msys.ElementForAbbreviation(a), i)
            self.assertEqual(msys.AbbreviationForElement(i), n)
            self.assertTrue(
                abs(msys.ElectronegativityForElement(i) - e) < 1e-6,
                msys.ElectronegativityForElement(i),
            )

    def testTopoIds(self):
        mol = msys.LoadDMS("tests/files/ww.dms")
        ids = msys.ComputeTopologicalIds(mol)
        self.assertEqual(len(ids), mol.natoms)

    def testGuessHydrogenPositions(self):
        mol = msys.LoadDMS("tests/files/ww.dms")
        hs = mol.select("hydrogen and not water")
        assert len(hs) > 10
        msys.GuessHydrogenPositions(hs)
        parents = [h.bonded_atoms[0] for h in hs]
        for h, p in zip(hs, parents):
            delta = h.pos - p.pos
            dist = (delta * delta).sum() ** 0.5
            self.assertTrue(dist > 0.8)
            self.assertTrue(dist < 1.2)

    def testAssignBondOrdersAndFormalCharges(self):
        # Smoke test only
        sys = msys.LoadDMS("tests/files/ww.dms")
        msys.AssignBondOrderAndFormalCharge(sys)
        msys.AssignBondOrderAndFormalCharge(sys.select("water"))

    def testAssignBondOrderTimeout(self):
        mol = msys.FromSmilesString(
            "c12c3c4c5c1c6c7c8c2c9c1c3c2c3c4c4c%10c5c5c6c6c7c7c%11c8c9c8c9c1c2c1c2c3c4c3c4c%10c5c5c6c6c7c7c%11c8c8c9c1c1c2c3c2c4c5c6c3c7c8c1c23"
        )
        with self.assertRaises(RuntimeError):
            msys.AssignBondOrderAndFormalCharge(mol, timeout=1)

        mol = msys.FromSmilesString("c1ccccc1")
        for b in mol.bonds:
            b.order = 1
        msys.AssignBondOrderAndFormalCharge(mol, timeout=1)
        orders = [b.order for b in mol.bonds]
        self.assertEqual(orders, [1, 2, 1, 2, 1, 2, 1, 1, 1, 1, 1, 1])

        for b in mol.bonds:
            b.order = 1
        msys.AssignBondOrderAndFormalCharge(mol)
        orders = [b.order for b in mol.bonds]
        self.assertEqual(orders, [1, 2, 1, 2, 1, 2, 1, 1, 1, 1, 1, 1])

    def testResonantCharges(self):
        mol = msys.Load("tests/files/jandor.sdf")
        msys.AssignBondOrderAndFormalCharge(mol)
        self.assertFalse("resonant_charge" in mol.atom_props)

        msys.AssignBondOrderAndFormalCharge(mol, compute_resonant_charges=True)
        topids = msys.ComputeTopologicalIds(mol)
        self.assertTrue("resonant_charge" in mol.atom_props)
        atms = [a for a in mol.atoms if a["resonant_charge"] != 0]
        aids = [topids[atm.id] for atm in atms]
        # print([a['resonant_charge'] for a in atms])
        # print([a.formal_charge for a in atms])
        # print(aids)
        # make sure we got a non-trivial resonant charge assignment
        assert aids[1] == aids[2]
        assert atms[1]["resonant_charge"] == atms[2]["resonant_charge"]
        assert atms[1].formal_charge != atms[2].formal_charge

        # exercise the other code paths
        msys.AssignBondOrderAndFormalCharge(mol.atoms, compute_resonant_charges=True)
        msys.AssignBondOrderAndFormalCharge(mol.atoms, 0, compute_resonant_charges=True)

    def testMemoryLeakPosVel(self):
        mol = msys.CreateSystem()
        for i in range(1000):
            mol.addAtom()

        pos = NP.zeros((mol.natoms, 3), "f")
        for i in range(100):
            mol.setPositions(pos)
            mol.setVelocities(pos)
        oldsize = vsize()
        for i in range(1000):
            mol.setPositions(pos)
            mol.setVelocities(pos)
        newsize = vsize()
        self.assertEqual(oldsize, newsize)

    def testBadPath(self):
        m = msys.CreateSystem()
        m.addAtom()
        with self.assertRaises(RuntimeError):
            msys.SaveDMS(m, "/root/nono.dms")
        with self.assertRaises(RuntimeError):
            msys.SaveMAE(m, "/root/nono.mae")
        with self.assertRaises(RuntimeError):
            msys.SavePDB(m, "/root/nono.pdb")

    def xxtestFailedWrite(self):
        print("building big system...")
        m = msys.CreateSystem()
        c = m.addChain()
        r = c.addResidue()
        id = r.id
        f = r._ptr.addAtom
        for i in range(1000000):
            f(id)
        print("saving...")
        msys.SaveMAE(m, "/usr/tmp/big.dms")
        print("done")

    def testAtom(self):

        m = msys.CreateSystem()
        a1 = m.addAtom()
        a2 = m.addAtom()
        self.assertTrue(a1.system == m)
        self.assertTrue(a2.system == m)
        self.assertEqual(a2.system.atoms[0], a1)
        self.assertEqual(a1.system.atoms[1], a2)
        self.assertNotEqual(m.atom(0), m.atom(1))

        self.assertNotEqual(a1, None)
        self.assertTrue(a1 != None)
        self.assertTrue(a1 == a1)
        self.assertFalse(a1 == None)
        self.assertFalse(a1 != a1)

        # FIXME: hmm, maybe fetching m.atom(0) ought to throw after the
        # atom has been removed.
        m.atom(0).remove()
        self.assertFalse(m.atom(0) in m.atoms)

        m2 = msys.CreateSystem()
        a3 = m2.addAtom()
        self.assertNotEqual(a1, a3)
        self.assertEqual(a3, a3)
        self.assertEqual(a3, m2.atom(0))

        res = a1.residue
        self.assertEqual(m.residues[0], res)
        rnew = m.addResidue()
        self.assertNotEqual(rnew, a1.residue)
        a1.residue = rnew
        self.assertEqual(rnew, a1.residue)

    def testSelectOnRemoved(self):
        m = msys.CreateSystem()
        m.addAtom()
        m.addAtom().atomic_number = 17
        m.addAtom()
        for a in m.select("atomicnumber 17"):
            a.remove()
        self.assertEqual(0, len(m.select("atomicnumber 17")))

    def testParamHandle(self):
        params = msys.CreateParamTable()
        self.assertNotEqual(params, None)
        self.assertEqual(params, params)
        self.assertNotEqual(params, msys.CreateParamTable())
        p1 = params.addParam()
        p2 = params.addParam()
        self.assertTrue(p1 == p1)
        self.assertTrue(p1 != p2)
        self.assertFalse(p1 != p1)
        self.assertFalse(p1 == p2)
        p1b = params.param(0)
        self.assertTrue(p1 == p1b)
        self.assertTrue(p1b != p2)
        self.assertFalse(p1b != p1)
        self.assertFalse(p1b == p2)

        self.assertEqual(len(set((p1, p2))), 2)
        self.assertEqual(len(set((p1, p1b))), 1)

    def testSystemHandle(self):
        m = msys.CreateSystem()
        self.assertNotEqual(m, None)
        self.assertEqual(m, m)
        self.assertNotEqual(m, msys.CreateSystem())

    def testTermHandle(self):
        m = msys.CreateSystem()
        a = m.addAtom()
        T = m.addTable("foo", 1)
        t1 = T.addTerm([a], None)
        t2 = T.addTerm([a], None)
        t1b = T.term(0)

        self.assertEqual(len(set((t1, t2))), 2)
        self.assertEqual(len(set((t1, t1b))), 1)

    def testCoalesce(self):
        m = msys.CreateSystem()
        a = m.addAtom()
        table = m.addTableFromSchema("posre_harm")
        p1 = table.params.addParam()
        p2 = table.params.addParam()
        p3 = table.params.addParam()
        p3["fcx"] = 32

        t1 = table.addTerm([a], p1)
        t2 = table.addTerm([a], p2)
        t3 = table.addTerm([a], p3)

        self.assertFalse(t1.param == t2.param)
        self.assertFalse(t1.param == t3.param)

        table.coalesce()
        self.assertTrue(t1.param == t2.param)
        self.assertTrue(t1.param == p1)
        self.assertTrue(t2.param == p1)
        self.assertFalse(t1.param == t3.param)
        self.assertTrue(t3.param == p3)

    def testAtomProps(self):
        m = msys.CreateSystem()
        a1 = m.addAtom()
        F, I, S = "F", "I", "S"
        for p, t in zip((F, I, S), (float, int, str)):
            m.addAtomProp(p, t)
            self.assertTrue(p in m.atom_props)
            self.assertEqual(m.atomPropType(p), t)
        self.assertEqual(len(m.atom_props), 3)
        a2 = m.addAtom()

        a1.charge = 3.5
        self.assertEqual(a1.charge, 3.5)
        self.assertEqual(m.atom(a1.id).charge, a1.charge)
        self.assertEqual(a1[F], 0)
        self.assertEqual(a1[I], 0)
        self.assertEqual(a1[S], "")

        a1[F] = 32.5
        a1[I] = 42
        a1[S] = "justinrocks"
        self.assertEqual(a1[F], 32.5)
        self.assertEqual(a1[I], 42)
        self.assertEqual(a1[S], "justinrocks")

        with self.assertRaises(KeyError):
            a1["foobar"] = 32

        self.assertTrue(hasattr(a2, "charge"))
        self.assertTrue("F" in a2)
        self.assertFalse("foobar" in a2)

        m.delAtomProp("F")
        self.assertFalse("F" in a2)
        self.assertEqual([m.atomPropType(n) for n in m.atom_props], [int, str])

    def testBond(self):
        m = msys.CreateSystem()
        a1 = m.addAtom()
        a2 = m.addAtom()
        a3 = m.addAtom()
        b = a2.addBond(a3)

        self.assertTrue(a1 != b)
        self.assertFalse(a1 == b)

        self.assertEqual(b.system, m)
        self.assertEqual(b.first, a2)
        self.assertEqual(b.second, a3)
        self.assertEqual(b.other(b.first), b.second)
        self.assertEqual(b.other(b.second), b.first)
        b.order = 32
        self.assertEqual(b.order, 32)
        self.assertEqual(len(m.bonds), 1)

        first, second = b.atoms
        self.assertEqual(first, a2)
        self.assertEqual(second, a3)

        self.assertEqual([a for a in a1.bonded_atoms], [])
        self.assertEqual([a for a in a2.bonded_atoms], [a3])
        self.assertEqual([a for a in a3.bonded_atoms], [a2])

        b.remove()
        self.assertEqual(len(m.bonds), 0)

        b1 = a1.addBond(a2)
        self.assertEqual(len(a2.bonds), 1)
        self.assertEqual(b1, m.findBond(a1, a2))
        self.assertEqual(b1, a1.findBond(a2))
        self.assertEqual(b1, a2.findBond(a1))
        self.assertEqual(None, m.findBond(a1, a3))

        self.assertEqual(a2.nbonds, 1)

        b.remove()
        b.remove()
        b.remove()
        self.assertEqual(len(a2.bonds), 1)

        with self.assertRaises(TypeError):
            m.delAtoms(m.bonds)
        with self.assertRaises(TypeError):
            m.delBonds(m.atoms)
        with self.assertRaises(TypeError):
            m.delResidues(m.chains)
        with self.assertRaises(TypeError):
            m.delChains(m.residues)

        m2 = msys.CreateSystem()
        with self.assertRaises(ValueError):
            m2.delAtoms(m.atoms)
        with self.assertRaises(ValueError):
            m2.delBonds(m.bonds)
        with self.assertRaises(ValueError):
            m2.delResidues(m.residues)
        with self.assertRaises(ValueError):
            m2.delChains(m.chains)

        m.delBonds(m.bonds)
        self.assertEqual(len(m.bonds), 0)
        m.delAtoms(m.atoms)
        self.assertEqual(len(m.atoms), 0)

        self.assertEqual(len(m.residues), 3)
        m.delResidues(m.residues[2:])
        self.assertEqual(len(m.residues), 2)
        m.delResidues(m.residues[0:1])
        self.assertEqual(len(m.residues), 1)
        m.delChains(m.chains)
        self.assertEqual(len(m.chains), 0)

    def testExportMAEContents(self):
        m = msys.Load("tests/files/noFused1.mae")
        tmp = tempfile.NamedTemporaryFile(suffix=".mae")
        path = tmp.name
        msys.Save(m, path)
        with open(path, "rb") as fp:
            b1 = fp.read()
        b2 = msys.SerializeMAE(m)
        self.assertEqual(b1, b2)

    def testExportDmsGz(self):
        old = msys.Load("tests/files/noFused1.mae")
        with tempfile.NamedTemporaryFile(suffix=".dms.gz") as tmp:
            msys.Save(old, tmp.name)
            new = msys.Load(tmp.name)
        self.assertEqual([a.name for a in old.atoms], [a.name for a in new.atoms])

    def testExportMaeGz(self):
        m = msys.Load("tests/files/noFused1.mae")
        tmp = tempfile.NamedTemporaryFile(suffix=".maegz")
        path = tmp.name
        f = gzip.open(path, "w")
        f.write(msys.SerializeMAE(m))
        f.write(msys.SerializeMAE(m))
        f.close()
        m2 = msys.Load(path)
        self.assertEqual(m2.ncts, 2)
        self.assertEqual(m2.natoms, 2 * m.natoms)

    def testRemove(self):
        m = msys.CreateSystem()
        r = m.addResidue()
        for i in range(10000):
            r.addAtom()
        self.assertEqual(r.natoms, 10000)
        r.remove()
        r.remove()
        r.remove()
        self.assertEqual(len(m.atoms), 0)

        c = m.addChain()
        for i in range(10000):
            c.addResidue()
        self.assertEqual(c.nresidues, 10000)

        c.remove()
        c.remove()
        c.remove()
        self.assertEqual(len(m.residues), 0)
        self.assertEqual(c.nresidues, 0)

    def testResidue(self):
        m = msys.CreateSystem()
        r = m.addResidue()
        a1 = r.addAtom()
        a2 = r.addAtom()
        self.assertEqual(r.system, a1.system)
        self.assertEqual(r.system, m)
        self.assertEqual(m, a2.system)
        self.assertEqual(len(m.residues), 1)
        self.assertEqual(len(m.atoms), 2)
        self.assertEqual(len(r.atoms), 2)
        r.name = "alA"
        r.resid = 32
        self.assertEqual(r.name, "alA")
        self.assertEqual(r.resid, 32)

        r.remove()
        r.remove()
        r.remove()
        self.assertEqual(len(m.residues), 0)
        self.assertEqual(len(m.atoms), 0)
        self.assertEqual(len(r.atoms), 0)

    def testSegid(self):
        m = msys.CreateSystem()
        with self.assertRaises(RuntimeError):
            m.addAtomProp("segid", str)
        r = m.addResidue()
        a1 = r.addAtom()
        a2 = r.addAtom()
        r.chain.segid = "WAT1"
        tmp = tempfile.NamedTemporaryFile(suffix=".dms")
        path = tmp.name
        msys.SaveDMS(m, path)
        m2 = msys.LoadDMS(path)
        msys.Load(path)
        self.assertEqual(1, m.nresidues)
        self.assertEqual(1, m2.nresidues)

    def testChain(self):
        m = msys.CreateSystem()
        c = m.addChain()
        r1 = c.addResidue()
        r2 = c.addResidue()
        a1 = r1.addAtom()
        a2 = r2.addAtom()
        self.assertEqual(c.system, r1.system)
        self.assertEqual(c.system, m)
        self.assertEqual(m, r2.system)
        self.assertEqual(len(m.chains), 1)
        self.assertEqual(len(m.residues), 2)
        self.assertEqual(len(c.residues), 2)
        self.assertEqual(len(m.atoms), 2)
        c.remove()
        self.assertEqual(len(m.residues), 0)
        self.assertEqual(len(c.residues), 0)
        self.assertEqual(len(m.atoms), 0)

    def testParamKeywords(self):
        m = msys.CreateSystem()
        T = m.addTable("foo", 1)
        P = T.params
        with self.assertRaises(KeyError):
            P.addParam(x=32)
        P.addProp("i", int)
        P.addProp("f", float)
        P.addProp("s", str)
        p = P.addParam(i=32, f=3.5, s="jrg")
        self.assertEqual(p["i"], 32)
        self.assertEqual(p["f"], 3.5)
        self.assertEqual(p["s"], "jrg")

    def testResetParams(self):
        m = msys.CreateSystem()
        a1 = m.addAtom()
        T = m.addTable("foo", 1)
        P = T.params
        P.addProp("f", float)
        p = P.addParam()
        p["f"] = 32
        t = T.addTerm([a1], p)

        P2 = msys.CreateParamTable()
        P2.addProp("i", int)
        p2 = P2.addParam()
        p2["i"] = 4

        with self.assertRaises(RuntimeError):
            T.term(0).param = p2

        T.params = P2
        self.assertEqual(T.params, P2)
        self.assertEqual(T.term(0).param, None)
        T.term(0).param = p2
        with self.assertRaises(RuntimeError):
            T.term(0).param = p

        T.params = P
        self.assertEqual(T.params, P)
        T.term(0).param = p

    def testTermTable(self):
        m = msys.CreateSystem()
        a1 = m.addAtom()
        a2 = m.addAtom()
        a3 = m.addAtom()
        a4 = m.addAtom()
        with self.assertRaises(RuntimeError):
            m.addTable("foo", 0)

        angle = m.addTable("angle", 3)

        # addTable with same name returns the same table
        angle2 = m.addTable("angle", 3)
        self.assertEqual(angle, angle2)

        # addTable with an explicit ParamTable must not already exist
        with self.assertRaises(RuntimeError):
            m.addTable("angle", 3, msys.CreateParamTable())

        self.assertEqual(angle.natoms, 3)
        self.assertEqual(len(angle.terms), 0)
        t1 = angle.addTerm(m.atoms[:3])

        self.assertTrue(angle.hasTerm(0))
        self.assertFalse(angle.hasTerm(1))

        self.assertEqual(t1.atoms, m.atoms[:3])
        self.assertTrue(t1.param is None)
        t2 = angle.addTerm(m.atoms[1:4], None)
        self.assertTrue(t2.param is None)
        self.assertEqual(len(angle.terms), 2)

        angle.name = "angle2"
        self.assertEqual(angle.name, "angle2")
        self.assertTrue("angle2" in m.table_names)
        self.assertFalse("angle" in m.table_names)
        angle.name = "angle"
        self.assertEqual(angle.name, "angle")
        self.assertTrue("angle" in m.table_names)
        self.assertFalse("angle2" in m.table_names)

        params = angle.params
        self.assertEqual(params, angle.params)
        p0 = params.addParam()
        p1 = params.addParam()
        self.assertEqual(p1.table, params)
        self.assertEqual(p0, params.param(0))
        self.assertEqual(p1, params.param(1))
        self.assertEqual(params.nparams, 2)
        t1.param = p1
        self.assertEqual(t1.param, p1)
        t1.param = None
        self.assertEqual(t1.param, None)

        with self.assertRaises(RuntimeError):
            angle.addTerm((m.atom(2),))
        with self.assertRaises(RuntimeError):
            angle.addTerm((m.atom(1), m.atom(1), m.atom(5)))

        angle.remove()
        self.assertFalse("angle" in m.table_names)
        self.assertFalse(angle in m.tables)

    def testIndex(self):
        m = msys.CreateSystem()
        a1 = m.addAtom()
        a2 = m.addAtom()
        a3 = m.addAtom()
        a4 = m.addAtom()
        t = m.addTable("foo", 2)

        self.assertEqual(t.findWithAll([a1]), [])
        self.assertEqual(t.findWithAny([a1]), [])

        t.addTerm((a1, a2))
        t.addTerm((a1, a3))
        t.addTerm((a1, a4))
        t.addTerm((a2, a3))
        t.addTerm((a4, a2))
        t.addTerm((a4, a3))
        self.assertEqual(t.nterms, 6)
        self.assertEqual(t.findWithAll([a1]), [t.term(x) for x in (0, 1, 2)])
        self.assertEqual(t.findWithAll([a2]), [t.term(x) for x in (0, 3, 4)])
        self.assertEqual(t.findWithAll([a3]), [t.term(x) for x in (1, 3, 5)])
        self.assertEqual(t.findWithAll([a4]), [t.term(x) for x in (2, 4, 5)])
        self.assertEqual(t.findWithAny([a1, a2]), [t.term(x) for x in range(5)])
        self.assertEqual(t.findWithAny([a1, a3]), [t.term(x) for x in (0, 1, 2, 3, 5)])

        self.assertEqual(t.findWithOnly([a1, a2, a4]), [t.term(x) for x in (0, 2, 4)])
        self.assertEqual(t.findWithOnly([a2, a1, a4]), [t.term(x) for x in (0, 2, 4)])
        self.assertEqual(
            t.findWithOnly([a4, a1, a2, a2]), [t.term(x) for x in (0, 2, 4)]
        )
        self.assertEqual(t.findWithOnly([]), [])
        self.assertEqual(t.findWithOnly([a4]), [])
        self.assertEqual(t.findWithOnly([a4]), [])

        self.assertEqual(t.findWithAll([a4, a2]), [t.term(4)])
        self.assertEqual(t.findWithAll([a2, a4]), [t.term(4)])
        self.assertEqual(t.findExact([a4, a2]), [t.term(4)])
        self.assertEqual(t.findExact([a2, a4]), [])
        a1.remove()

        self.assertEqual(t.findWithOnly([a2, a4]), [t.term(4)])
        self.assertEqual(t.findWithAll([a1]), [t.term(x) for x in ()])
        self.assertEqual(t.findWithAll([a2]), [t.term(x) for x in (3, 4)])
        self.assertEqual(t.findWithAll([a3]), [t.term(x) for x in (3, 5)])
        self.assertEqual(t.findWithAll([a4]), [t.term(x) for x in (4, 5)])
        a5 = m.addAtom()
        t.addTerm((a5, a3))
        self.assertEqual(t.findWithAll([a5]), [t.term(6)])
        self.assertEqual(t.findWithAll([a3]), [t.term(x) for x in (3, 5, 6)])
        self.assertEqual(t.findWithAny([a3, a5]), [t.term(x) for x in (3, 5, 6)])
        t.addTerm((a5, a5))
        self.assertEqual(t.findWithAny([a5]), [t.term(x) for x in (6, 7)])
        self.assertEqual(t.findWithAll([a5]), [t.term(x) for x in (6, 7)])

    def testParamIndex(self):
        params = msys.CreateParamTable()
        params.addProp("i", int)
        params.addProp("f", float)
        params.addProp("s", str)
        p0 = params.addParam()
        p1 = params.addParam()
        p2 = params.addParam()
        p3 = params.addParam()
        p0["i"] = 1
        p1["i"] = 1
        p1["f"] = 1.5
        p2["f"] = 1.5
        p2["s"] = "x"
        p3["s"] = "x"
        for (k, v), ids in {
            ("i", 0): (2, 3),
            ("i", 1): (0, 1),
            ("f", 0): (0, 3),
            ("f", 1.5): (1, 2),
            ("s", ""): (0, 1),
            ("s", "x"): (2, 3),
        }.items():
            self.assertEqual(params.find(k, v), [params.param(x) for x in ids])
        for i in range(100):
            params.addParam()
        for (k, v), ids in {
            ("i", 1): (0, 1),
            ("f", 1.5): (1, 2),
            ("s", "x"): (2, 3),
        }.items():
            self.assertEqual(params.find(k, v), [params.param(x) for x in ids])
        p3["s"] = "y"
        for (k, v), ids in {
            ("i", 1): (0, 1),
            ("f", 1.5): (1, 2),
            ("s", "x"): (2,),
        }.items():
            self.assertEqual(params.find(k, v), [params.param(x) for x in ids])

    def testDelTerm(self):
        m = msys.CreateSystem()
        a1 = m.addAtom()
        a2 = m.addAtom()
        a3 = m.addAtom()
        a4 = m.addAtom()
        table = m.addTable("foo", 2)
        table.addTerm((a1, a2))
        table.addTerm((a1, a3))
        table.addTerm((a1, a4))
        table.addTerm((a2, a3))
        table.addTerm((a4, a2))
        table.addTerm((a4, a3))
        self.assertEqual(table.nterms, 6)
        table.term(2).remove()  # removes term 2
        self.assertEqual(table.nterms, 5)

        table.delTermsWithAtom(a3)  # removes term 1,3,5
        self.assertEqual(table.nterms, 2)
        a4.remove()
        self.assertEqual(table.nterms, 1)

    def testParamProps(self):
        params = msys.CreateParamTable()
        self.assertEqual(len(params.props), 0)
        for n, t in zip(("F", "I", "S"), (float, int, str)):
            params.addProp(n, t)
            self.assertTrue(n in params.props)
            self.assertEqual(params.propType(n), t)
        self.assertEqual(params.props, ["F", "I", "S"])

        params.delProp("I")
        self.assertEqual(params.props, ["F", "S"])

        p1 = params.addParam()
        p1["F"] = 1.5
        p2 = p1.duplicate()
        self.assertEqual(p1["F"], p2["F"])
        p2["F"] = 2.5
        self.assertNotEqual(p1["F"], p2["F"])

    def testDirectParamProps(self):
        """just like testParamProps, but operate directly on the TermTable
        instead of its ParamTable"""
        m = msys.CreateSystem()
        table = m.addTable("table", 5)
        self.assertEqual(len(table.term_props), 0)
        for n, t in zip(("F", "I", "S"), (float, int, str)):
            table.addTermProp(n, t)
            self.assertTrue(n in table.term_props)
            self.assertEqual(table.termPropType(n), t)

    def testPropOverlap(self):
        """ ensure term props cannot overlap param props. """
        m = msys.CreateSystem()
        table = m.addTable("table", 1)
        table.params.addProp("a", float)
        with self.assertRaises(RuntimeError):
            table.addTermProp("a", int)

    def testParamWithNoProps(self):
        m = msys.CreateSystem()
        a = m.addAtom()
        table = m.addTable("foo", 1)
        p = table.params.addParam()
        table.addTerm([a], p)
        table.category = "bond"
        self.assertEqual(table.category, "bond")
        tmp = tempfile.NamedTemporaryFile(suffix=".dms")
        path = tmp.name
        msys.SaveDMS(m, path)

    def testSaveDmsStructureOnly(self):
        m = msys.CreateSystem()
        a = m.addAtom()
        table = m.addTable("foo", 1)
        p = table.params.addParam()
        table.addTerm([a], p)
        table.category = "bond"
        self.assertEqual(table.category, "bond")
        tmp = tempfile.NamedTemporaryFile(suffix=".dms")
        path = tmp.name
        msys.SaveDMS(m, path, structure_only=False)
        self.assertEqual(len(msys.Load(path).tables), 1)
        msys.Save(m, path, structure_only=False)
        self.assertEqual(len(msys.Load(path).tables), 1)
        msys.SaveDMS(m, path, structure_only=True)
        self.assertEqual(len(msys.Load(path).tables), 0)
        msys.Save(m, path, structure_only=True)
        self.assertEqual(len(msys.Load(path).tables), 0)

    def testLoadStructureOnly(self):
        m = msys.CreateSystem()
        m.addAtom()
        m.addAtom().atomic_number = 6
        for ext in ".dms", ".mae":
            with tempfile.NamedTemporaryFile(suffix=ext) as fp:
                path = fp.name
                msys.Save(m, path)
                mol = msys.Load(path)
                self.assertEqual(mol.natoms, 2, "w/ structure for %s" % ext)
                mol = msys.Load(path, structure_only=True)
                self.assertEqual(mol.natoms, 1, "w/o structure for %s" % ext)

                plugin = molfile.guess_filetype(path)
                mol = plugin.read(path)
                self.assertEqual(mol.natoms, 2, "molfile for %s" % ext)

    def testLoadWithoutTables(self):
        mol = msys.CreateSystem()
        a0 = mol.addAtom()
        a1 = mol.addAtom()
        a2 = mol.addAtom()
        a3 = mol.addAtom()

        a0.atomic_number = 6
        a1.atomic_number = 6
        a2.atomic_number = 6
        a3.atomic_number = 0

        a0.addBond(a1)
        a0.addBond(a2)

        t = mol.addTableFromSchema("stretch_harm")
        p = t.params.addParam()
        t.addTerm([a1, a2], p)

        # TODO add support for prmtop (use tests/files/sys.prmtop as input)
        for suffix in ".dms", ".mae":
            with tempfile.NamedTemporaryFile(suffix=suffix) as fp:
                path = fp.name
                msys.Save(mol, path)
                m = msys.Load(path)
                self.assertEqual(m.natoms, 4)
                self.assertEqual(m.nbonds, 2)
                t = m.table("stretch_harm")
                self.assertEqual(t.nterms, 1)
                self.assertEqual([a.id for a in t.term(0).atoms], [1, 2])

                m = msys.Load(path, without_tables=True)
                self.assertEqual(m.natoms, 4)
                self.assertEqual(m.nbonds, 2)
                t = m.getTable("stretch_harm")
                self.assertEqual(t, None)

                m = msys.Load(path, structure_only=True, without_tables=False)
                self.assertEqual(m.natoms, 3)
                self.assertEqual(m.nbonds, 2)
                t = m.getTable("stretch_harm")
                self.assertEqual(t.nterms, 1)
                self.assertEqual([a.id for a in t.term(0).atoms], [1, 2])

                m = msys.Load(path, structure_only=True)
                self.assertEqual(m.natoms, 3)
                self.assertEqual(m.nbonds, 2)
                t = m.getTable("stretch_harm")
                self.assertEqual(t, None)

    def testFunnyNames(self):
        tmp = tempfile.NamedTemporaryFile(suffix=".dms")
        path = tmp.name
        m = msys.CreateSystem()
        aux = msys.CreateParamTable()
        aux.addProp("group", float)
        aux.addParam()
        m.addAuxTable("aux", aux)
        msys.SaveDMS(m, path)

        m2 = m.clone()
        aux2 = m2.auxtable("aux")
        self.assertEqual(aux.nprops, aux2.nprops)
        self.assertEqual(aux.nparams, aux2.nparams)

        m3 = msys.CreateSystem()
        m3.append(m)
        aux3 = m3.auxtable("aux")
        self.assertEqual(aux.nprops, aux3.nprops)
        self.assertEqual(aux.nparams, aux3.nparams)

        table = m.addTable("foo", 1)
        table.category = "bond"
        table.addTermProp("select", int)
        table.params.addProp("group", float)
        m.addAtom()
        table.addTerm(m.atoms, table.params.addParam())
        msys.SaveDMS(m, path)

    def testFormalCharge(self):
        tmp = tempfile.NamedTemporaryFile(suffix=".dms")
        path = tmp.name
        m = msys.CreateSystem()
        a = m.addAtom()
        self.assertEqual(a.formal_charge, 0)
        a.formal_charge = 32
        self.assertEqual(a.formal_charge, 32)
        a.formal_charge = -10
        self.assertEqual(a.formal_charge, -10)
        msys.SaveDMS(m, path)
        m2 = msys.LoadDMS(path)
        self.assertEqual(m.atom(0).formal_charge, -10)

    def testLoadDmsBuffer(self):
        with open("tests/files/2f4k.dms", "rb") as fp:
            s = fp.read()
        msys.LoadDMS(buffer=s)

    def testLoadDmsMany(self):
        tmp = tempfile.NamedTemporaryFile(suffix=".dms")
        for mol in (
            msys.Load("tests/files/%s" % m)
            for m in ("1vcc.mae", "cofactors.sdf", "ww.dms")
        ):
            tmp.write(msys.FormatDMS(mol))
        tmp.flush()
        mols = [m for m in msys.LoadMany(tmp.name)]
        self.assertEqual(len(mols), 3)
        self.assertEqual([m.natoms for m in mols], [26228, 562, 7614])

    def testRefcount(self):
        m = msys.CreateSystem()
        a = m.addAtom()
        table = m.addTable("foo", 1)
        ptr = table.params._ptr
        p1 = table.params.addParam()
        p2 = table.params.addParam()
        p3 = table.params.addParam()
        t1 = table.addTerm([a], p1)
        self.assertEqual(ptr.refcount(t1.param.id), 1)
        t2 = table.addTerm([a], p1)
        self.assertEqual(ptr.refcount(t1.param.id), 2)
        self.assertEqual(ptr.refcount(t2.param.id), 2)
        t1.remove()
        self.assertEqual(ptr.refcount(t2.param.id), 1)
        t3 = table.addTerm([a], p2)
        self.assertEqual(ptr.refcount(t2.param.id), 1)
        self.assertEqual(ptr.refcount(t3.param.id), 1)

        t3.remove()
        self.assertEqual(ptr.refcount(t2.param.id), 1)

    def testSharedParams(self):
        m1 = msys.CreateSystem()
        m2 = msys.CreateSystem()
        m1.addAtom()
        m2.addAtom()
        m2.addAtom()
        params = msys.CreateParamTable()
        p1 = params.addParam()
        p2 = params.addParam()
        table1 = m1.addTable("table", 1, params)
        table2 = m2.addTable("table", 1, params)
        self.assertEqual(table1.params, table2.params)
        t1 = table1.addTerm(m1.atoms, p2)
        t2 = table2.addTerm(m2.atoms[1:], p2)
        self.assertEqual(t1.param, t2.param)
        self.assertEqual(t2.param, p2)

        params.addProp("fc", float)
        p1["fc"] = 32
        p2["fc"] = 42
        self.assertEqual(t1.param["fc"], 42)
        self.assertEqual(t2.param["fc"], 42)

        # refcounts work even across shared tables
        t1["fc"] = 52
        self.assertEqual(t2["fc"], 42)

    def testTermParamProps(self):
        m = msys.CreateSystem()
        m.addAtom()
        table = m.addTable("table", 1)
        p = table.params.addParam()
        t = table.addTerm(m.atoms, p)
        t2 = table.addTerm(m.atoms, p)

        with self.assertRaises(KeyError):
            t["x"]
            t.stateB["x"]
            t.stateA["x"]

        table.params.addProp("x", float)
        table.params.addProp("y", str)
        self.assertEqual(t["x"], 0)
        self.assertEqual(t["y"], "")

        t2["x"] = 32
        t2["y"] = "whodat"
        self.assertEqual(t["x"], 0)
        self.assertEqual(t["y"], "")
        self.assertEqual(table.params.nparams, 2)

        self.assertEqual(t2._ptr.getProp(t2.id, 0), 32)
        t2._ptr.setProp(t2.id, 0, 33)
        self.assertEqual(t2._ptr.getProp(t2.id, 0), 33)
        t2._ptr.setProp(t2.id, 0, 32)

    def testTermProps(self):
        m = msys.CreateSystem()
        m.addAtom()
        table = m.addTable("table", 1)
        F = "F"
        I = "I"
        table.addTermProp(F, float)
        t1 = table.addTerm(m.atoms)
        table.addTermProp(I, int)
        t2 = table.addTerm(m.atoms)
        self.assertEqual(t1[F], 0)
        self.assertEqual(t1[I], 0)
        self.assertEqual(t2[F], 0)
        self.assertEqual(t2[I], 0)
        t1[F] = 3.5
        t2[I] = 42
        self.assertEqual(t1[F], 3.5)
        self.assertEqual(t1[I], 0)
        self.assertEqual(t2[F], 0)
        self.assertEqual(t2[I], 42)
        self.assertEqual(table.term_props, [F, I])

        m2 = msys.CreateSystem()
        m2.addAtom()
        with self.assertRaises(RuntimeError):
            table.addTerm(m2.atoms)
        alist = m.atoms
        alist[0] = m2.atoms[0]
        with self.assertRaises(RuntimeError):
            table.addTerm(alist)
        alist[0] = m.atoms[0]
        table.addTerm(alist)

        table.delTermProp(I)
        with self.assertRaises(KeyError):
            t1[I]

    def testAux(self):
        m = msys.CreateSystem()
        e = msys.CreateParamTable()
        m.addAuxTable("cmap", e)
        self.assertTrue("cmap" in m.auxtable_names)
        self.assertEqual(e, m.auxtable("cmap"))
        self.assertTrue(e in m.auxtables)
        m.delAuxTable("cmap")
        self.assertFalse("cmap" in m.auxtable_names)
        self.assertFalse(e in m.auxtables)

    def testSchemas(self):
        tmp = tempfile.NamedTemporaryFile(suffix=".dms")
        path = tmp.name
        m = msys.CreateSystem()
        for s in msys.TableSchemas():
            m.addTableFromSchema(s)
        for s in msys.NonbondedSchemas():
            m.addNonbondedFromSchema(s)
            m.nonbonded_info.vdw_funct = ""
        msys.Save(m, path)
        msys.Load(path)

    def testSchemaWithName(self):
        m = msys.CreateSystem()
        t = m.addTableFromSchema("dihedral_trig", "my_name")

    def testSavePartialNonbonded(self):
        m = msys.CreateSystem()
        a = m.addAtom()
        t = m.addNonbondedFromSchema("vdw_12_6", "arithemetic/geometric")
        with tmpfile(suffix=".dms") as tmp:
            # don't allow partial nonbonded
            with self.assertRaises(RuntimeError):
                msys.Save(m, tmp.name)
            # ... unless it's structure-only
            msys.Save(m, tmp.name, structure_only=True)
            # fixed the nonbonded table
            t.addTerm([a], t.params.addParam())
            msys.Save(m, tmp.name, structure_only=False)
            msys.Save(m, tmp.name, structure_only=True)

    def testTableDMS(self):
        assert "stretch_harm" in msys.TableSchemas()
        assert "vdw_12_6" in msys.NonbondedSchemas()

        m = msys.CreateSystem()
        t = m.addTableFromSchema("stretch_harm")
        self.assertEqual(t.natoms, 2)
        self.assertEqual(t.params.props, ["r0", "fc"])
        self.assertEqual([t.params.propType(n) for n in t.params.props], [float, float])
        self.assertEqual(t.term_props, ["constrained"])
        self.assertEqual(t.termPropType("constrained"), int)

    def testNonbondedDMS(self):
        m = msys.CreateSystem()
        nb = m.addNonbondedFromSchema("vdw_12_6", "arithemetic/geometric")
        nb2 = m.addNonbondedFromSchema("vdw_12_6", "arithemetic/geometric")
        m.addNonbondedFromSchema("vdw_12_6")
        self.assertEqual(nb, nb2)
        params = nb.params
        props = params.props
        self.assertEqual(props, ["sigma", "epsilon"])
        self.assertEqual([nb.params.propType(n) for n in props], [float, float])

    def testMultipleNonbonded(self):
        m = msys.CreateSystem()
        a = m.addAtom()
        m.addAtom()
        disp = m.addTable("nonbonded_dispersion", 1)
        repl = m.addTable("nonbonded_repulsion", 1)
        elec = m.addTable("nonbonded_charge", 1)
        for t in disp, repl, elec:
            t.category = "nonbonded"
            p = t.params.addParam()
            t.addTerm([a], p)

        disp.params.addProp("foo", float)
        repl.params.addProp("bar", float)
        elec.params.addProp("charge", float)
        m.nonbonded_info.vdw_funct = "disp_repl_charge"
        m.nonbonded_info.vdw_rule = "geom/geom/geom"

        tmp = tempfile.NamedTemporaryFile(suffix=".dms")
        path = tmp.name

        try:
            msys.SaveDMS(m, path)
        finally:
            pass
        t = m.addTable("nonbonded", 1)
        t.category = "nonbonded"
        p = t.params.addParam()
        t.addTerm([a], p)
        t.addTerm([m.atom(1)], p)
        try:
            msys.SaveDMS(m, path)
        finally:
            pass

    def testGlobalCell(self):
        m = msys.CreateSystem()
        m.cell[0][0] = 32
        self.assertEqual(m.cell[0][0], 32)
        tgt = [1, 2, 3]
        m.cell[0][:] = tgt
        m.cell[1][:] = tgt
        self.assertEqual(m.cell[0].tolist(), m.cell[1].tolist())
        self.assertNotEqual(m.cell[0].tolist(), m.cell[2].tolist())
        with self.assertRaises(IndexError):
            m.cell[3]
        with self.assertRaises(IndexError):
            m.cell[-4]

        c = m.getCell()
        self.assertEqual(c.tolist(), m.cell.tolist())
        c[1][2] = 99
        self.assertNotEqual(c.tolist(), m.cell.tolist())
        m.setCell(c)
        self.assertEqual(c.tolist(), m.cell.tolist())
        m.setCell(c.astype("f"))
        self.assertEqual(c.tolist(), m.cell.tolist())

        # check lifetime
        m2 = msys.CreateSystem()
        c = m2.cell
        del m2
        c
        c
        self.assertEqual(c.tolist(), [[0, 0, 0], [0, 0, 0], [0, 0, 0]])

    def testNonbondedInfo(self):
        m = msys.CreateSystem()
        nb = m.nonbonded_info
        nb.vdw_funct = "justinrocks"
        nb.vdw_rule = "yep"
        nb.es_funct = "oh yeah"
        self.assertEqual(nb.vdw_funct, "justinrocks")
        self.assertEqual(nb.vdw_rule, "yep")
        self.assertEqual(nb.es_funct, "oh yeah")

    def testAppendNonbondedInfo(self):
        A = msys.CreateSystem()
        B = msys.CreateSystem()
        attrs = dict(vdw_funct="vf", vdw_rule="vr", es_funct="es")
        for k, v in attrs.items():
            setattr(B.nonbonded_info, k, v)

        A.append(B)
        for k, v in attrs.items():
            self.assertEqual(getattr(B.nonbonded_info, k), v)
            self.assertEqual(getattr(A.nonbonded_info, k), v)

    def testCloneShareParams(self):
        m = msys.CreateSystem()
        a1 = m.addAtom()
        t1 = m.addTable("t1", 1)
        t1.category = "nonbonded"
        t1.params.addParam()
        t1.addTerm([a1], t1.params.param(0))

        m2 = m.clone()
        self.assertNotEqual(t1.params, m2.table("t1").params)
        m2 = m.clone(share_params=True)
        self.assertEqual(t1.params, m2.table("t1").params)

    def testCloneNumpyIntegers(self):
        m = msys.CreateSystem()
        m.addAtom()
        m.clone(NP.arange(1, dtype="i"))

    def testClone2(self):
        m = msys.CreateSystem()
        a1 = m.addAtom()
        a2 = m.addAtom()
        t1 = m.addTable("t1", 1)
        t1.category = "nonbonded"
        params = t1.params
        t2 = m.addTable("t2", 1, params)
        t2.category = "nonbonded"

        p1 = params.addParam()
        p2 = params.addParam()
        t1.addTerm([a1], p1)
        t1.addTerm([a2], p2)
        t2.addTerm([a2], p1)

        self.assertEqual(params.nparams, 2)

        m2 = m.clone()
        t1 = m2.table("t1")
        t2 = m2.table("t2")
        params = t1.params
        self.assertEqual(t1.params, t2.params)
        self.assertEqual(params.nparams, 2)

    def testClone(self):
        m = msys.CreateSystem()
        m.addAtom().name = "a"
        m.addAtom().name = "b"
        m.addAtom().name = "c"
        b1 = m.atom(0).addBond(m.atom(1))
        with self.assertRaises(KeyError):
            b1["res_order"] = 3.5
        m.addBondProp("res_order", float)
        b1["res_order"] = 3.5

        c = m.clone(m.atoms[::-1])
        self.assertEqual([a.name for a in c.atoms], ["c", "b", "a"])
        self.assertEqual(c.bond(0)["res_order"], 3.5)

        with self.assertRaises(RuntimeError):
            m.clone([m.atoms[0], m.atoms[1], m.atoms[0]])

        m.atom(2).remove()
        with self.assertRaises(RuntimeError):
            m.clone([m.atom(0), m.atom(1), m.atom(2)])

    def testKeys(self):
        m = msys.CreateSystem()
        T = m.addTable("foo", 1)
        P = T.params
        P.addProp("x", float)
        T.addTermProp("y", int)

        a0 = m.addAtom()
        p = T.params.addParam()
        t = T.addTerm([a0], p)
        self.assertEqual(list(t.keys()), ["x", "y"])
        self.assertEqual(list(p.keys()), ["x"])

    def testAppend(self):
        m = msys.CreateSystem()
        a0 = m.addAtom()
        a0.name = "a"
        a1 = m.addAtom()
        a1.name = "b"
        b = a0.addBond(a1)
        m.addAtomProp("foo", float)
        m.addBondProp("bar", int)
        a1["foo"] = 3.14
        b["bar"] = 42
        m.append(m)
        self.assertEqual([a.name for a in m.atoms], ["a", "b", "a", "b"])
        self.assertEqual(m.atom(2)["foo"], 0)
        self.assertEqual(m.atom(3)["foo"], 3.14)
        self.assertEqual(m.bond(1)["bar"], 42)

    def testReplaceWithSortedTerms(self):
        mol = msys.Load("tests/files/1vcc.mae")
        perm = NP.random.permutation(mol.natoms).tolist()
        new = mol.clone(perm)
        self.assertFalse(
            all(t.id == t.atoms[0].id for t in new.table("nonbonded").terms)
        )
        params1 = {t.atoms[0]: t["sigma"] for t in new.table("nonbonded").terms}
        new.table("nonbonded").replaceWithSortedTerms()
        self.assertTrue(
            all(t.id == t.atoms[0].id for t in new.table("nonbonded").terms)
        )
        params2 = {t.atoms[0]: t["sigma"] for t in new.table("nonbonded").terms}
        self.assertEqual(params1, params2)

    def testRemoveTable(self):
        m = msys.CreateSystem()
        a = m.addAtom()
        T = m.addTable("foo", 1)
        T.params.addProp("x", int)
        p = T.params.addParam()
        t1 = T.addTerm([a], p)
        t1["x"] = 32
        self.assertEqual(t1.param.id, 0)

        T2 = m.addTable("bar", 1, T.params)
        t2 = T2.addTerm([a], p)
        t2["x"] = 32
        self.assertEqual(t2.param.id, 1)
        self.assertEqual(T.params.nparams, 2)

        T.name = "fooj"
        self.assertEqual(T.system, m)
        self.assertEqual(T.name, "fooj")
        T.remove()
        self.assertEqual(T.system, None)
        self.assertEqual(T.name, "")
        with self.assertRaises(RuntimeError):
            T.addTerm([a], p)
        self.assertEqual(T.name, "")
        self.assertEqual([t.id for t in T.terms], [])

    def testChargeC(self):
        m = msys.CreateSystem()
        nb = m.addNonbondedFromSchema("vdw_12_6")
        alc = m.addTable("alchemical_nonbonded", 1, nb.params)
        alc.category = "nonbonded"
        alc.addTermProp("chargeC", float)

        a = m.addAtom()
        p = nb.params.addParam()

        nb.addTerm([a], p)
        t = alc.addTerm([a], p)
        t["chargeC"] = 0.5

        tmp = tempfile.NamedTemporaryFile(suffix=".dms")
        path = tmp.name
        msys.SaveDMS(m, path)
        m2 = msys.LoadDMS(path)
        msys.SaveDMS(m2, path)
        alc = m2.table("alchemical_nonbonded")
        self.assertEqual(alc.term(0)["chargeC"], 0.5)

    def testUpdateFragids(self):
        m = msys.CreateSystem()
        a1 = m.addAtom()
        a2 = m.addAtom()
        a3 = m.addAtom()

        def fragids(m):
            return [a.fragid for a in m.atoms]

        frags = m.updateFragids()
        self.assertEqual(fragids(m), [0, 1, 2])
        self.assertEqual(frags, [[m.atom(i)] for i in range(3)])

        a1.addBond(a3)
        frags = m.updateFragids()
        self.assertEqual(fragids(m), [0, 1, 0])
        self.assertEqual(frags, [[m.atom(0), m.atom(2)], [m.atom(1)]])

    def testBadDMS(self):
        tmp = tempfile.NamedTemporaryFile(suffix=".dms")
        path = tmp.name
        m = msys.CreateSystem()
        msys.SaveDMS(m, path)
        a0 = m.addAtom()
        msys.SaveDMS(m, path)
        nb = m.addNonbondedFromSchema("vdw_12_6", "geometric")
        # can't save - no terms for particle 0
        with self.assertRaises(RuntimeError):
            msys.SaveDMS(m, path)
        t = nb.addTerm([a0])
        # can't save - no param
        with self.assertRaises(RuntimeError):
            msys.SaveDMS(m, path)

        p0 = nb.params.addParam()
        t.param = p0
        # now we can save
        msys.SaveDMS(m, path)
        # twiddle param to something bad
        t._ptr.setParam(t.id, None)
        with self.assertRaises(RuntimeError):
            msys.SaveDMS(m, path)
        t._ptr.setParam(t.id, 0)

        stretch = m.addTableFromSchema("stretch_harm")
        t = stretch.addTerm([a0, a0])
        with self.assertRaises(RuntimeError):
            msys.SaveDMS(m, path)
        p = stretch.params.addParam()
        t.param = p
        msys.SaveDMS(m, path)
        t._ptr.setParam(t.id, None)
        with self.assertRaises(RuntimeError):
            msys.SaveDMS(m, path)

        # the failed saves deleted the tmpfile, which causes a warning
        # message when the tmpfile object is closed.  Suppress that.
        try:
            tmp.close()
        except:
            pass

    def testPositions(self):
        m = msys.CreateSystem()
        a0 = m.addAtom()
        a1 = m.addAtom()
        a2 = m.addAtom()
        a0.y = 1
        a1.z = 2
        a2.x = 3
        p = m.positions
        self.assertEqual(p[0][1], 1)
        self.assertEqual(p[1][2], 2)
        self.assertEqual(p[2][0], 3)
        p[1][2] = 4
        m.positions = p
        p = m.positions
        self.assertEqual(p[1][2], 4)
        self.assertTrue((p == m.positions).all())
        with self.assertRaises(ValueError):
            m.positions = p[1:, :]
        with self.assertRaises(ValueError):
            m.positions = p[:, 1:]

        # test of the low level Id-based accessors for positions
        ids = m._ptr.selectAsList("index 1")
        m._ptr.setPositions([[8, 9, 10]], ids)
        self.assertEqual(list(m.positions[1]), [8, 9, 10])
        self.assertEqual(list(m._ptr.getPositions(ids)[0]), [8, 9, 10])

    def testVelocities(self):
        m = msys.CreateSystem()
        a0 = m.addAtom()
        a1 = m.addAtom()
        a2 = m.addAtom()
        a0.vy = 1
        a1.vz = 2
        a2.vx = 3
        p = m.getVelocities()
        self.assertEqual(p[0][1], 1)
        self.assertEqual(p[1][2], 2)
        self.assertEqual(p[2][0], 3)
        p[1][2] = 4
        m.setVelocities(p)
        p = m.getVelocities()
        self.assertEqual(p[1][2], 4)
        self.assertTrue((p == m.getVelocities()).all())
        with self.assertRaises(ValueError):
            m.setVelocities(p[1:, :])
        with self.assertRaises(ValueError):
            m.setVelocities(p[:, 1:])

        # test of the low level Id-based accessors for velocities
        ids = m._ptr.selectAsList("index 1")
        m._ptr.setVelocities([[8, 9, 10]], ids)
        self.assertEqual(list(m.getVelocities()[1]), [8, 9, 10])
        self.assertEqual(list(m._ptr.getVelocities(ids)[0]), [8, 9, 10])

    def testSelectChain(self):
        m = msys.CreateSystem()
        c0 = m.addChain()
        c1 = m.addChain()
        c2 = m.addChain()
        c0.name = "A"
        c1.name = "A"
        c1.segid = "PRO"
        c2.name = "B"
        self.assertEqual(m.selectChain("B"), c2)
        self.assertEqual(m.selectChain(segid="PRO"), c1)
        self.assertEqual(m.selectChain("A", "PRO"), c1)
        self.assertEqual(m.selectChain("C"), None)
        self.assertEqual(m.selectChain(segid="foo"), None)
        with self.assertRaises(ValueError):
            m.selectChain("A")
        c0.remove()
        c2.remove()
        self.assertEqual(m.selectChain(), c1)
        c1.remove()
        self.assertEqual(m.selectChain(), None)

    def testSelectResidue(self):
        m = msys.CreateSystem()
        c = m.addChain()
        r0 = c.addResidue()
        r1 = c.addResidue()
        r2 = c.addResidue()
        r0.name = "A"
        r1.name = "A"
        r1.resid = 1
        r2.name = "B"
        r2.insertion = "A"
        self.assertEqual(c.selectResidue(name="B"), r2)
        self.assertEqual(c.selectResidue(resid=1), r1)
        self.assertEqual(c.selectResidue(1, "A"), r1)
        self.assertEqual(c.selectResidue(3), None)
        self.assertEqual(c.selectResidue(resid=4), None)
        with self.assertRaises(ValueError):
            c.selectResidue(name="A")
        r0.remove()
        r2.remove()
        self.assertEqual(c.selectResidue(), r1)
        r1.remove()
        self.assertEqual(c.selectResidue(), None)

    def testSelectAtom(self):
        m = msys.CreateSystem()
        r = m.addResidue()
        a0 = r.addAtom()
        a1 = r.addAtom()
        a2 = r.addAtom()
        a0.name = "C"
        a1.name = "CB"
        a2.name = "C"
        self.assertEqual(r.selectAtom("CB"), a1)
        self.assertEqual(r.selectAtom("CA"), None)
        with self.assertRaises(ValueError):
            r.selectAtom("C")

    def testSelectCt(self):
        m = msys.CreateSystem()
        c0 = m.addCt()
        c1 = m.addCt()
        c2 = m.addCt()
        c0.name = "C"
        c1.name = "CB"
        c2.name = "C"
        self.assertEqual(m.selectCt("CB"), c1)
        self.assertEqual(m.selectCt("CA"), None)
        with self.assertRaises(ValueError):
            m.selectCt("C")

    def testConcatenatedMae(self):
        d = os.path.dirname(__file__)
        m = msys.Load(os.path.join(d, "files/two.mae"))
        self.assertEqual(m.ncts, 4)
        self.assertEqual(m.nchains, 4)
        self.assertEqual(m.natoms, 9498)

    def testMaeWhitespaceKey(self):
        m = msys.CreateSystem()
        m.addAtom()
        m.ct(0)["key with whitespace"] = 42
        msys.SerializeMAE(m)

    def testAppendMae(self):
        tmp = tempfile.NamedTemporaryFile(suffix=".mae")
        path = tmp.name
        m = msys.CreateSystem()
        m.addAtom().atomic_number = 1

        msys.SaveMAE(m, path, append=True)
        msys.SaveMAE(m, path, append=True)
        m2 = msys.Load(path)
        self.assertEqual(m2.natoms, 2)

        tmp = tempfile.NamedTemporaryFile(suffix=".mae")
        path = tmp.name
        msys.SaveMAE([m, m], path, append=True)
        m2 = msys.Load(path)
        self.assertEqual(m2.natoms, 2)

    def testAppendMol2(self):
        m = msys.CreateSystem()
        m.addAtom().atomic_number = 1

        tmp = tempfile.NamedTemporaryFile(suffix=".mol2")
        path = tmp.name
        m.ct(0).name = "XYZ"
        msys.SaveMol2(m, path, append=True)
        m.ct(0).name = "ABC"
        msys.SaveMol2(m, path, append=True)
        for i, m in enumerate(msys.LoadMany(path)):
            self.assertEqual(m.natoms, 1)
            self.assertEqual(m.ct(0).name, ("XYZ", "ABC")[i])
            self.assertEqual(m.name, ("XYZ", "ABC")[i])
        self.assertEqual(i, 1)

    def testMol2Residues(self):
        m1 = msys.Load("tests/files/test_UID_not_reordered.mol2")
        m2 = msys.Load("tests/files/test_UID_corrected.mol2")
        self.assertEqual(m1.nresidues, m2.nresidues)

    def testMol2Arginine(self):
        arg = msys.Load("tests/files/arginine.pdb")
        msys.AssignBondOrderAndFormalCharge(arg)
        tmp = tempfile.NamedTemporaryFile(suffix=".mol2")
        path = tmp.name
        msys.SaveMol2(arg, path)

        with open(path) as fp:
            data = fp.read()
        assert " N.pl3 " in data, data
        assert " ar\n" in data, data

    def testAppendSDF(self):
        m = msys.CreateSystem()
        m.addAtom().atomic_number = 1
        tmp = tempfile.NamedTemporaryFile(suffix=".sdf")

        msys.Save(m, tmp.name, append=True)
        msys.Save(m, tmp.name, append=True)
        for i, m in enumerate(msys.LoadMany(tmp.name)):
            self.assertEqual(m.natoms, 1)
        self.assertEqual(i, 1)

    def testAlchemicalMaeRestraint(self):
        d = os.path.dirname(__file__)
        m = msys.Load(os.path.join(d, "files/alchemical_restraint.mae"))
        self.assertEqual(m.natoms, 5)
        con = m.table("constraint_ah2")
        res = m.table("posre_harm")
        self.assertEqual(con.nterms, 1)
        self.assertEqual(sorted(x.id for x in con.term(0).atoms), [0, 1, 2])
        self.assertEqual(res.nterms, 2)
        self.assertEqual(res.term(0).atoms[0].id, 0)
        self.assertEqual(res.term(1).atoms[0].id, 1)
        self.assertEqual(res.term(0)["fcx"], 0.25)
        self.assertEqual(res.term(1)["fcy"], 0.35)

    def testMaeFBHW(self):
        d = os.path.dirname(__file__)
        m = msys.Load(os.path.join(d, "files/fbhw.mae"))
        posre = m.table("posre_fbhw")
        self.assertEqual(posre.params.nparams, 1)
        self.assertEqual(posre.nterms, 32)

    def testRingSystems(self):
        for path, match_systems in {
            "tests/files/fused.sdf": [[0, 1]],
            "tests/files/noFused1.mae": [[0]],
            "tests/files/noFused2.mae": [[0], [1], [2]],
            "tests/files/jandor.sdf": [[0], [1], [2, 3]],
            "tests/files/colzuy.mae.gz": [[0, 1, 2]],
            "tests/files/kanzoo.mae.gz": [[0, 1], [2, 3]],
        }.items():

            mol = msys.Load(path)
            systems = msys.GetRingSystems(mol.atoms)
            self.assertEqual(sorted(systems), match_systems)

    def testSSSR(self):
        sys = msys.Load("tests/files/cubane.dms")
        rings = msys.GetSSSR(sys.atoms, True)
        self.assertTrue(len(rings) == 6)
        for ring in rings:
            self.assertTrue(len(ring) == 4)
            for ring2 in rings:
                if ring != ring2:
                    intersect = set([a.id for a in ring]) & set([a.id for a in ring2])
                    self.assertTrue(len(intersect) == 2 or len(intersect) == 0)
        r0 = [r for r in rings if sys.atom(0) in r]
        self.assertTrue(len(r0) == 3)

    def testAnnotatedSystem(self):
        # Test aromaticity and atom props
        sys = msys.CreateSystem()
        res = sys.addChain().addResidue()
        c = [res.addAtom() for i in range(6)]
        h = [res.addAtom() for i in range(6)]
        cc = []
        ch = []
        for i in range(6):
            c[i].atomic_number = 6
            cc.append(c[i].addBond(c[(i + 1) % 6]))
            h[i].atomic_number = 1
            ch.append(h[i].addBond(c[i]))
        with self.assertRaises(RuntimeError):
            annot_sys = msys.AnnotatedSystem(sys)
        msys.AssignBondOrderAndFormalCharge(sys)
        annot_sys = msys.AnnotatedSystem(sys)
        repr(annot_sys)
        for i in range(6):
            self.assertTrue(annot_sys.aromatic(c[i]))
            self.assertTrue(annot_sys.aromatic(cc[i]))
            self.assertTrue(not annot_sys.aromatic(h[i]))
            self.assertTrue(not annot_sys.aromatic(ch[i]))
            self.assertTrue(annot_sys.hcount(c[i]) == 1)
            self.assertTrue(annot_sys.degree(c[i]) == 3)
            self.assertTrue(annot_sys.valence(c[i]) == 4)
            self.assertTrue(annot_sys.loneelectrons(c[i]) == 0)
            self.assertTrue(annot_sys.hybridization(c[i]) == 2)
            self.assertTrue(annot_sys.ringbondcount(c[i]) == 2)
            self.assertTrue(annot_sys.hcount(h[i]) == 0)
            self.assertTrue(annot_sys.degree(h[i]) == 1)
            self.assertTrue(annot_sys.valence(h[i]) == 1)
            self.assertTrue(annot_sys.loneelectrons(h[i]) == 0)
            self.assertTrue(annot_sys.ringbondcount(h[i]) == 0)
        h[5].remove()
        o = res.addAtom()
        o.atomic_number = 8
        co = c[5].addBond(o)
        h[5] = res.addAtom()
        h[5].atomic_number = 1
        oh = h[5].addBond(o)
        msys.AssignBondOrderAndFormalCharge(sys)
        annot_sys = msys.AnnotatedSystem(sys)
        for i in range(6):
            self.assertTrue(annot_sys.aromatic(c[i]))
            self.assertTrue(annot_sys.aromatic(cc[i]))
        self.assertTrue(not annot_sys.aromatic(co))
        self.assertTrue(not annot_sys.aromatic(o))
        new_h = [res.addAtom() for i in range(2)]
        for i in range(2):
            new_h[i].atomic_number = 1
            new_h[i].addBond(c[i])
        msys.AssignBondOrderAndFormalCharge(sys)
        annot_sys = msys.AnnotatedSystem(sys)
        for i in range(6):
            self.assertTrue(not annot_sys.aromatic(c[i]))
            self.assertTrue(not annot_sys.aromatic(cc[i]))

    def testHashAnnotatedSystem(self):
        m1 = msys.CreateSystem()
        m1.addAtom().atomic_number = 6
        a1 = msys.AnnotatedSystem(m1)

        m2 = msys.CreateSystem()
        m2.addAtom().atomic_number = 6
        a2 = msys.AnnotatedSystem(m2)
        self.assertTrue(a1 == a1)
        self.assertTrue(a2 == a2)
        self.assertTrue(a2 != a1)
        self.assertTrue(a1 != a2)
        self.assertTrue(hash(a1) != hash(a2))
        self.assertTrue(len(set((a1, a2))) == 2)

    def testSmartsPattern(self):
        import ast

        d = os.path.dirname(__file__)
        files = [
            "ww.dms",
            "membrane.dms",
            "4A9C_chainA_46.mae",
            "4A9C_chainA_cf_21.mae",
            "acrd.mae",
            "azul.mae",
            "boron.mae",
            "fenz.mae",
            "indz.mae",
            "ndph.mae",
            "pegm.mae",
            "PC2777373.mae",
        ]
        for f in files:
            name = f[:-4]
            mol = msys.Load("tests/smarts_tests/" + f, structure_only=True)
            msys.AssignBondOrderAndFormalCharge(mol)
            annot_mol = msys.AnnotatedSystem(mol)
            atoms = mol.select("not water")
            with open("tests/smarts_tests/%s_matches" % name) as fp:
                tests = ast.literal_eval(fp.read())
            for k, v in tests.items():
                sp = msys.SmartsPattern(k)
                match = sp.findMatches(annot_mol, atoms)
                if match != v:
                    msg = (
                        "SMARTS mismatch for molecule '%s', pattern '%s'.\nMissing"
                        " matches: %s\nExtra matches: %s\n"
                    )
                    msg = msg % (
                        name,
                        k,
                        str([i for i in v if i not in match]),
                        str([i for i in match if i not in v]),
                    )
                    self.assertTrue(False, msg)

    def testGraphColors(self):
        G = msys.Graph
        mol1 = msys.Load("tests/files/tip5p.mae").clone("fragid 0")
        mol2 = mol1.clone()

        # 1) no color: pseudos ignored, heavy atoms matched 1-1
        m1 = G(mol1).match(G(mol2))
        self.assertEqual(sorted(k.id for k in m1), [0, 1, 2])
        for k, v in m1.items():
            self.assertEqual(k.id, v.id)

        # 2) color from atomic number: same result
        anum = lambda x: x.atomic_number
        m2 = G(mol1, anum).match(G(mol2, anum))
        self.assertEqual(m1, m2)

        # 3) color from atomic number + 1: match everything 1-1 since pseudos are no longer ignored
        anum = lambda x: x.atomic_number + 1
        m3 = G(mol1, anum).match(G(mol2, anum))
        self.assertEqual(sorted(k.id for k in m3), [0, 1, 2, 3, 4])
        for k, v in m3.items():
            self.assertEqual(k.id, v.id)

        # 4) use particle name as the attribute for vsites; vsites are Vrt0 and Vrt1
        from collections import defaultdict

        d = defaultdict(lambda: len(d))
        colors = [
            (a.atomic_number if a.atomic_number > 0 else d[a.name] + 128)
            for a in mol1.atoms
        ]
        self.assertEqual(colors, [8, 1, 1, 128, 129])
        m4 = G(mol1, colors).matchAll(G(mol2, colors))
        self.assertEqual(len(m4), 2)
        colors[-1] = colors[-2]
        m5 = G(mol1, colors).matchAll(G(mol2, colors))
        self.assertEqual(len(m5), 4)

    def testGraphDistinctAtoms(self):
        mol = msys.CreateSystem()
        mol.addAtom().atomic_number = 2
        with self.assertRaises(RuntimeError):
            msys.Graph([mol.atom(0), mol.atom(0)])

    def testGraph(self):
        sys_A = msys.CreateSystem()
        res = sys_A.addResidue()
        a0 = res.addAtom()
        a1 = res.addAtom()
        a2 = res.addAtom()
        a3 = res.addAtom()
        a4 = res.addAtom()
        a5 = res.addAtom()
        a6 = res.addAtom()
        a7 = res.addAtom()
        a8 = res.addAtom()
        a0.addBond(a5)
        a0.addBond(a6)
        a0.addBond(a7)
        a0.addBond(a8)
        a1.addBond(a2)
        a2.addBond(a3)
        a2.addBond(a4)
        a3.addBond(a6)
        a4.addBond(a5)
        a0.atomic_number = 4
        a1.atomic_number = 1
        a2.atomic_number = 4
        a3.atomic_number = 3
        a4.atomic_number = 2
        a5.atomic_number = 4
        a6.atomic_number = 4
        a7.atomic_number = 1
        a8.atomic_number = 1
        atoms_A = [a0, a2, a3, a4, a5, a6, a7, a8]
        sys_B = msys.CreateSystem()
        res = sys_B.addResidue()
        b0 = res.addAtom()
        b1 = res.addAtom()
        b2 = res.addAtom()
        b3 = res.addAtom()
        b4 = res.addAtom()
        b5 = res.addAtom()
        b6 = res.addAtom()
        b7 = res.addAtom()
        b8 = res.addAtom()
        b0.addBond(b1)
        b0.addBond(b2)
        b1.addBond(b3)
        b1.addBond(b5)
        b1.addBond(b7)
        b2.addBond(b6)
        b4.addBond(b6)
        b4.addBond(b7)
        b6.addBond(b8)
        b0.atomic_number = 4
        b1.atomic_number = 4
        b2.atomic_number = 3
        b3.atomic_number = 1
        b4.atomic_number = 2
        b5.atomic_number = 1
        b6.atomic_number = 4
        b7.atomic_number = 4
        b8.atomic_number = 1
        atoms_B = [b0, b1, b2, b3, b4, b5, b6, b7]
        graph_A = msys.Graph(atoms_A)
        graph_B = msys.Graph(atoms_B)
        self.assertEqual(sorted(atoms_B), sorted(graph_B.atoms()))
        self.assertTrue(graph_A != graph_B)
        # Test match
        matches = graph_A.match(graph_B)
        match1 = {a0: b1, a2: b6, a3: b2, a4: b4, a5: b7, a6: b0, a7: b3, a8: b5}
        match2 = {a0: b1, a2: b6, a3: b2, a4: b4, a5: b7, a6: b0, a7: b5, a8: b3}
        self.assertTrue(matches == match1 or matches == match2)
        all_matches = graph_A.matchAll(graph_B)
        self.assertTrue(len(all_matches) == 2)
        self.assertTrue(match1 in all_matches and match2 in all_matches)
        a0.atomic_number = 3
        graph_A = msys.Graph(atoms_A)
        # Test no match because of atomic number
        self.assertTrue(graph_A.match(graph_B) is None)
        self.assertTrue(len(graph_A.matchAll(graph_B)) == 0)
        a0.atomic_number = 4
        bond = sys_A.findBond(a1, a2)
        sys_A.delBonds([bond])
        graph_A = msys.Graph(atoms_A)
        # Test no match because of bond topology
        self.assertTrue(graph_A.match(graph_B) is None)
        self.assertTrue(len(graph_A.matchAll(graph_B)) == 0)

        # Test matchAll on cube graph---8 atoms connected as a cube, plus
        # two extra atoms attached to each of these 8 vertices
        sys_A = msys.CreateSystem()
        res = sys_A.addResidue()
        a = [res.addAtom() for i in range(24)]
        for i in range(24):
            a[i].atomic_number = 1
        for i in range(4):
            a[i].addBond(a[(i + 1) % 4])
            a[i + 4].addBond(a[(i + 1) % 4 + 4])
        for i in range(4):
            a[i].addBond(a[i + 4])
        for i in range(8):
            a[i].addBond(a[8 + 2 * i])
            a[i].addBond(a[8 + 2 * i + 1])
        graph_A = msys.Graph(sys_A.atoms)
        all_matches = graph_A.matchAll(graph_A)
        self.assertTrue(len(all_matches) == 8 * 6 * 2 ** 8)

        # Test matchAll substructure search
        graph_B = msys.Graph([a[0], a[1], a[3], a[4]])
        all_matches = graph_B.matchAll(graph_A, substructure=True)
        self.assertTrue(len(all_matches) == 8 * 6)
        self.assertTrue(len(all_matches[0]) == 4)


class TestWrap(unittest.TestCase):
    def setUp(self):
        bonds = (
            (0, 1),
            (0, 2),
            (4, 3),
            (4, 6),
            (5, 6),
            (7, 8),
            (9, 8),
            (9, 10),
            (10, 11),
            (11, 7),
        )

        agg = (1, 8, 9)

        box = NP.array(
            [[0.70710678, 0.70710678, 0.0], [-1.41421356, 1.41421356, 0.0], [0, 0, 3]]
        )

        pos = NP.array(
            [
                # first fragment
                0,
                2,
                0,
                0,
                4.1,
                0,
                0.1,
                2.2,
                0.3,
                # second fragment */
                3,
                1,
                3,
                3,
                1.1,
                3,
                3.1,
                1.2,
                3.3,
                3.2,
                1.0,
                6.5,
                # third fragment */
                1,
                0,
                2,
                1,
                0.8,
                2,
                1.1,
                4.2,
                2.3,
                0.8,
                0.0,
                2.5,
                0.9,
                0.0,
                2.1,
            ],
            "f",
        )

        pos = pos.reshape((12, 3))
        top = [[] for _ in range(len(pos))]
        for i, j in bonds:
            top[i].append(j)
            top[j].append(i)

        self.top = top
        self.box = box
        self.pos = pos
        self.agg = agg

    def testBondsGlueAggWrap(self):
        p = pfx.Pfx(self.top, fixbonds=True)
        p.glue(self.agg)
        p.apply(self.pos, self.box)
        NP.testing.assert_almost_equal(
            self.box,
            [
                [0.7071067690849304, 0.7071067690849304, 0.0],
                [-1.4142135381698608, 1.4142135381698608, 0.0],
                [0.0, 0.0, 3.0],
            ],
        )
        NP.testing.assert_almost_equal(
            self.pos,
            [
                [-1.4142135381698608, 0.5857864618301392, 0.0],
                [-0.7071067690849304, 0.564466118812561, 0.0],
                [-1.314213514328003, 0.785786509513855, 0.30000001192092896],
                [-0.5355339050292969, 0.29289329051971436, 0.0],
                [-0.5355339050292969, 0.39289331436157227, 0.0],
                [-0.4355340003967285, 0.4928933382034302, 0.2999999523162842],
                [-0.33553385734558105, 0.29289329051971436, 0.5],
                [0.2928932309150696, -0.7071067690849304, -1.0],
                [-0.41421353816986084, -0.6142135262489319, -1.0],
                [-0.31421345472335815, -0.042640864849090576, -0.7000000476837158],
                [0.09289324283599854, -0.7071067690849304, -0.5],
                [0.19289320707321167, -0.7071067690849304, -0.9000000953674316],
            ],
        )


class TestImporter(unittest.TestCase):
    def test_1vcc(self):
        mol = msys.Load("tests/files/1vcc.mae")
        pro = mol.clone("protein")
        new = mol.clone("protein and noh")
        imp = msys.SystemImporter(new)
        imp.initialize(new.atoms)
        for a in mol.select("protein and hydrogen"):
            imp.addAtom(
                a.residue.chain.name,
                a.residue.chain.segid,
                a.residue.resid,
                a.residue.name,
                a.name,
            )
        self.assertEqual(new.natoms, pro.natoms)
        self.assertEqual(new.nresidues, pro.nresidues)


class TestSvd(unittest.TestCase):
    def testNice(self):
        A = NP.array([[1, 2, 3], [3, 2, 0], [4, 1, 5],])

        u, w, vt = NP.linalg.svd(A)
        U, W, V = pfx.svd_3x3(A)
        u, w, vt = NP.linalg.svd(A)
        U, W, V = pfx.svd_3x3(A)
        u, w, vt = NP.linalg.svd(A)
        U, W, V = pfx.svd_3x3(A)
        t0 = time()
        u, w, vt = NP.linalg.svd(A)
        t1 = time()
        U, W, V = pfx.svd_3x3(A)
        t2 = time()
        print("\nNP %8.3fms PFX %8.3fms" % ((t1 - t0) * 1000, (t2 - t1) * 1000))
        u[:, 0] *= -1
        u[:, 2] *= -1
        v = vt.transpose()
        v[:, 0] *= -1
        v[:, 2] *= -1
        NP.testing.assert_almost_equal(w, W)
        NP.testing.assert_almost_equal(u, U)
        NP.testing.assert_almost_equal(v, V)

    def testSingularRows(self):
        A = NP.array([[1, 2, 3], [1, 2, 3], [4, 1, 5],])

        u, w, vt = NP.linalg.svd(A)
        U, W, V = pfx.svd_3x3(A)
        u[:, 0] *= -1
        v = vt.transpose()
        v[:, 0] *= -1
        v[:, 2] *= -1
        NP.testing.assert_almost_equal(w, W)
        NP.testing.assert_almost_equal(u, U)
        NP.testing.assert_almost_equal(v, V)

    def testSingularCols(self):
        A = NP.array([[1, 1, 3], [3, 3, 0], [4, 4, 5],])

        u, w, vt = NP.linalg.svd(A)
        U, W, V = pfx.svd_3x3(A)
        u[:, 0] *= -1
        v = vt.transpose()
        v[:, 0] *= -1
        NP.testing.assert_almost_equal(w, W)
        NP.testing.assert_almost_equal(u, U)
        NP.testing.assert_almost_equal(v, V)


if __name__ == "__main__":
    unittest.main(verbosity=2)
