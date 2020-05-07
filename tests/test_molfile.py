import os, sys, unittest

import sqlite3

import msys
from msys import molfile
from msys.molfile import findframe
import unittest
import numpy
import shutil as SH
import subprocess
import tempfile

# we import pands here because it maybe imported later during a
# test, which results in a weird, spurious error involving importlib.
import pandas

'''
FIXME: This is sort of a mess.  There are tests for commonly used plugins,
the generic python interface, and the caching infrastructure all jumbled
together here.  
'''

class TestXyz(unittest.TestCase):
    def test1(self):
        dms = molfile.dms.read('tests/files/ch4.dms')
        dtr = molfile.dtr.read('tests/files/ch4.dtr')
        with tempfile.NamedTemporaryFile(suffix='.xyz') as tmp:
            xyz = molfile.xyz.write(tmp.name, atoms=dms.atoms)
            xyz.frame(dms.frame(0))
            xyz.frame(dtr.frame(0))
            xyz.close()
            r = molfile.xyz.read(tmp.name)
            self.assertEqual(r.natoms, 5)
            self.assertEqual(r.nframes, -1)
            n=0
            for f in r.frames():
                n+=1
            self.assertEqual(n,2)


class TestDcd(unittest.TestCase):
    def setUp(self):
        self.r = molfile.dcd.read("tests/files/alanin.dcd")

    def testNframes(self):
        self.assertEqual(self.r.nframes, 100)

    def testNatoms(self):
        self.assertEqual(self.r.natoms, 66)

    def testRandomAccess(self):
        import random
        random.seed(1492)
        fids=list(range(100))
        random.shuffle(fids)
        frames=[f for f in self.r.frames()]
        for fid in fids:
            self.assertTrue((self.r.frame(fid).pos == frames[fid].pos).all())

class GuessFiletypeTestCase(unittest.TestCase):
    def testDtr(self):
        p=molfile.guess_filetype("foo.dtr")
        self.assertEqual(p.prettyname, "DESRES Trajectory, clobber")

class DtrTestCase(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.PATH = os.path.join(tempfile.mkdtemp(), 'bad.dtr')

    def setUp(self):
        ''' creates a dtr with metadata, a timekeys file with just a header,
        and not_hashed/.ddparams '''
        SH.rmtree(self.PATH, ignore_errors=True)
        # get the dtr writer to use 1 frame per file
        self.natoms=10000000
        self.writer=molfile.dtr.write(self.PATH, natoms=self.natoms)

    def testReadEmptyDtr(self):
        self.writer.close()
        assert molfile.DtrReader(self.PATH).nframes == 0

    def testCheckpointBytes(self):
        kv = molfile.DtrReader("/f/a/jobstep/14388415/0/checkpoint.atr").keyvals(0)
        self.assertEqual(kv["TITLE"], "ANTON")
        self.assertEqual(kv["ANTON2_CHECKPOINT_CM_TIMES_SQUARE"][:16], b'\x1d\x00\x00\x00Times square')

        CHK = 'ANTON2_CHECKPOINT_ITEMS'
        blob = kv[CHK]
        tmpdir = tempfile.mkdtemp(dir='/tmp', suffix='.dtr')
        writer = msys.molfile.DtrWriter(tmpdir, natoms=10)
        writer.append(kv['CHEMICALTIME'][0], kv)
        writer.sync()
        
        atr1 = molfile.DtrReader(tmpdir)
        blob1 = atr1.keyvals(0)[CHK]
        
        self.assertEqual(blob, blob1)

    def testDtrBox(self):
        for dbl in False, True:
            self.setUp()
            f=molfile.Frame(self.natoms, double_precision=dbl)
            for i in range(3): f.box[i][i]=i+1
            self.writer.frame(f)
            self.writer.sync()
            r=molfile.dtr.read(self.PATH)
            for f in r.frames(): 
                for i in range(3):
                    self.assertEqual(f.box[i][i], i+1)

    def tearDown(self):
        del self.writer
        SH.rmtree(self.PATH, ignore_errors=True)

    def addFrame(self, time):
        f=molfile.Frame(self.natoms, False)
        f.time=time
        self.writer.frame(f)
        self.writer.sync()

    def testEmpty(self):
        r=molfile.DtrReader(self.PATH)
        self.assertEqual(0,r.natoms)
        self.assertEqual(0,r.nframes)
        with self.assertRaises(IOError): r.frame(0)

        dtr=molfile.dtr.read(self.PATH)
        self.assertEqual(dtr.at_time_near(32), None)
        self.assertEqual(dtr.at_time_lt(32), None)
        self.assertEqual(dtr.at_time_le(32), None)
        self.assertEqual(dtr.at_time_gt(32), None)
        self.assertEqual(dtr.at_time_ge(32), None)

    def testListFields(self):
        k1 = molfile.list_fields('tests/files/force_groups.dtr')
        self.assertEqual(k1[0], 'dim__alchemical_angle_harm')
        self.assertEqual(len(k1), 45)

        k2 = molfile.list_fields('tests/files/ch4.dtr')
        self.assertEqual(k2, ['CHEMICAL_TIME', 'ENERGY', 'EX_ENERGY', 'FORMAT', 'KIN_ENERGY', 'POSITION', 'POT_ENERGY', 'PRESSURE', 'PRESSURETENSOR', 'TEMPERATURE', 'TITLE', 'UNITCELL', 'VIRIALTENSOR'])
        self.assertEqual(len(k2), 13)

        k3 = molfile.list_fields('tests/files/relative.stk')
        self.assertEqual(k3, k2)

        etr = '/d/vault/dhmvault-12/anton2/125/9576573/0000/energy.etr'
        if os.path.exists(etr):
            k4 = molfile.list_fields(etr)
            self.assertTrue('TSS_UH' in k4)
            self.assertEqual(len(k4), 21)
        else:
            assert False, "No energy.etr"


    def testBadFrame(self):
        self.addFrame(1.0)
        self.addFrame(3.0)
        r=molfile.DtrReader(self.PATH)
        r.frame(0)
        with self.assertRaises(IOError): r.frame(2)
        with self.assertRaises(IOError): r.frame(3)

    def testBadFrame2(self):
        self.addFrame(1.0)
        r=molfile.DtrReader(self.PATH)
        r.frame(0)
        with self.assertRaises(IOError): r.frame(1)

    def testNoTimekeys(self):
        os.unlink('%s/timekeys' % self.PATH)
        with self.assertRaises(RuntimeError): molfile.DtrReader(self.PATH)

    def testShortTimekeys(self):
        open('%s/timekeys' % self.PATH, 'w').close()
        with self.assertRaises(RuntimeError): molfile.DtrReader(self.PATH)

    def testBadTimekeys(self):
        with open('%s/timekeys' % self.PATH, 'w') as fd: fd.write('hello')
        with self.assertRaises(RuntimeError): molfile.DtrReader(self.PATH)

    def testNoMeta(self):
        ''' test missing metadata frame.  DESRESCode#1525 '''
        # Need at least one frame in the dtr or metadata gets ignored 
        self.addFrame(0.0)
        os.unlink('%s/metadata' % self.PATH)
        with self.assertRaises(RuntimeError):
            molfile.DtrReader(self.PATH)

    def testOneFramePerFile(self):
        ''' test reading dtr with one frame per file.  DESRESCode#1525 '''
        self.addFrame(0.0)
        molfile.DtrReader(self.PATH).frame(0)
        self.addFrame(1.0)
        molfile.DtrReader(self.PATH).frame(0)

    def testMetadataFields(self):
        self.addFrame(0.0)
        meta = molfile.DtrReader(self.PATH).metadata
        self.assertEqual(meta['WORKDIR'], os.getcwd())
        self.assertEqual(meta['MSYS_VERSION'], msys.version.version)


    def testTruncatedFrame(self):
        self.addFrame(0.0)
        self.addFrame(1.0)
        self.addFrame(2.0)
        r=molfile.DtrReader(self.PATH)
        open(r.fileinfo(1)[0], 'w').close()
        r.frame(0)
        with self.assertRaises(IOError): r.frame(1)
        r.frame(2)

    def testZeroedFrame(self):
        ''' fill a frame file with zeros '''
        self.addFrame(0.0)
        self.addFrame(1.0)
        self.addFrame(2.0)
        r=molfile.DtrReader(self.PATH)
        info = r.fileinfo(1)
        path = info[0]
        sz=info[3]
        with open(path, 'w') as fd: print('\0' * sz, file=fd)
        r.frame(0)
        with self.assertRaises(IOError): r.frame(1)
        r.frame(2)

    def testFrameKeywords(self):
        ''' construct Frame in various ways '''
        with self.assertRaises(Exception):
            f=molfile.Frame()       # natoms is required
        f=molfile.Frame(natoms=32)
        self.assertEqual(len(f.pos), 32)
        self.assertEqual(f.vel, None)
        self.assertEqual(f.dpos, None)
        self.assertEqual(f.dvel, None)

        f=molfile.Frame(natoms=32, double_precision=True)
        self.assertEqual(len(f.dpos), 32)
        self.assertEqual(f.vel, None)
        self.assertEqual(f.pos, None)
        self.assertEqual(f.dvel, None)

        f=molfile.Frame(natoms=32, with_velocities=True)
        self.assertEqual(len(f.pos), 32)
        self.assertEqual(len(f.vel), 32)
        self.assertEqual(f.dpos, None)
        self.assertEqual(f.dvel, None)

        f=molfile.Frame(natoms=32, with_velocities=True, double_precision=True)
        self.assertEqual(len(f.dpos), 32)
        self.assertEqual(len(f.dvel), 32)
        self.assertEqual(f.pos, None)
        self.assertEqual(f.vel, None)

class TestStk(unittest.TestCase):
    STK='tests/files/run.stk'

    def testEmptyFramesetAtStartOfStk(self):
        s = molfile.DtrReader('/d/en/fep-6/2020/02/27/gullingj_k0kvwpri/pyanton/PA104_0/energy_globals.stk')
        self.assertEqual(s.nframes, 10412)

    def testMixedAtomCount(self):
        with self.assertRaises(RuntimeError):
            s1 = molfile.DtrReader('tests/files/mixed_atom_count1.stk')
        with self.assertRaises(RuntimeError):
            s1 = molfile.DtrReader('tests/files/mixed_atom_count2.stk')

    def testLsFile(self):
        s = molfile.StkFile.read('tests/files/ene.ls')
        self.assertEqual(s.nframes, 4167)

    def testRelativePath(self):
        molfile.dtr.read('tests/files/relative.stk').frame(0)

    def testReplaceWithDifferentAtoms(self):
        ''' create an stk, read from it to create a cache entry, then 
        overwrite it with a new one that has a different number of atoms.
        There should be no issues. '''

        #os.environ['DTRPLUGIN_VERBOSE']='1'
        tmp=tempfile.mkdtemp()
        m=molfile.dtr.write('%s/1.dtr' % tmp, natoms=10)
        m.frame(molfile.Frame(10, False))
        m.close()
        junkfile = tempfile.NamedTemporaryFile(suffix='.stk')
        junk = junkfile.name
        with open(junk, 'w') as f:
            print('%s/1.dtr' % tmp, file=f)

        #print "----- first read of junk.stk ------"
        r=molfile.dtr.read(junk)
        self.assertEqual(r.nframes, 1)
        self.assertEqual(r.natoms, 10)
        self.assertEqual(r.has_velocities, False)
        self.assertTrue(r.frame(0).vel is None)

        #print "----- second read of junk.stk ------"
        r=molfile.dtr.read(junk)
        self.assertEqual(r.nframes, 1)
        self.assertEqual(r.natoms, 10)
        self.assertEqual(r.has_velocities, False)
        self.assertTrue(r.frame(0).vel is None)

        #print "----- modify junk.stk ------"

        m=molfile.dtr.write('%s/2.dtr' % tmp, natoms=20)
        f = molfile.Frame(20, True)
        f.time=0
        m.frame(f)
        f.time=1
        m.frame(f)
        m.close()
        with open(junk, 'w') as f:
            print('%s/2.dtr' % tmp, file=f)

        #print "----- read modified junk.stk ------"
        r=molfile.dtr.read(junk)
        self.assertEqual(r.nframes, 2)
        self.assertEqual(r.natoms, 20)
        self.assertEqual(r.has_velocities, True)
        self.assertTrue(r.frame(0).vel is not None)

        ## remove the velocities again
        m=molfile.dtr.write('%s/3.dtr' % tmp, natoms=20)
        f = molfile.Frame(20, False)
        f.time=0
        m.frame(f)
        f.time=1
        m.frame(f)
        m.close()
        with open(junk, 'w') as f:
            print('%s/3.dtr' % tmp, file=f)
        r=molfile.dtr.read(junk)
        self.assertEqual(r.nframes, 2)
        self.assertEqual(r.natoms, 20)
        self.assertEqual(r.has_velocities, False)
        self.assertTrue(r.frame(0).vel is None)

        SH.rmtree(tmp)


    @unittest.skipIf(os.getenv('DESRES_LOCATION')!='EN', 'Runs only from EN location')
    def testTimes(self):
        stk=molfile.dtr.read(self.STK)
        times=stk.times
        for i in 0,1,100,1000,10000, 25000:
            self.assertEqual(times[i], stk.at_time_near(times[i]).time)

    @unittest.skipIf(os.getenv('DESRES_LOCATION')!='EN', 'Runs only from EN location')
    def testDtrReader(self):
        dtr=molfile.DtrReader(self.STK)
        sizemap=dict()
        for i in range(dtr.nframesets):
            sizemap[dtr.frameset_path(i)]=dtr.frameset_size(i)
            self.assertEqual( dtr.frameset_is_compact(i), True)

        countmap=dict()
        for i in range(dtr.nframes):
            dtrpath=dtr.fileinfo(i)[7]
            countmap.setdefault(dtrpath,0)
            countmap[dtrpath] += 1
        self.assertEqual( sizemap, countmap) 

        times=dtr.times()
        for i in 0,1,100,1000,10000, 25000:
            self.assertEqual(i, dtr.index_near(times[i]))

    @unittest.skipIf(os.getenv('DESRES_LOCATION')!='EN', 'Runs only from EN location')
    def testStk(self):

        frames = 1, 10, 100,1000, 10000
        dtr=molfile.DtrReader(self.STK)
        self.assertEqual(dtr.nframes, 27159)
        
        for i in frames:
          path, time, offset, framesize, first, last, filesize, dtrpath, dtrsize = dtr.fileinfo(i)
          with open(path, 'rb') as fp:
              allbytes = fp.read()
          self.assertEqual(len(allbytes), filesize)
          bytes=allbytes[offset:offset+framesize]
        
          frame=dtr.frame(i, bytes)
          p1 = frame.pos
          p2 = dtr.frame(i).pos
          self.assertTrue((p1==p2).all())
          v1 = frame.vel
          v2 = dtr.frame(i).vel
          self.assertTrue(v1 is not None)
          self.assertTrue(v2 is not None)
          self.assertTrue((v1==v2).all())
        
class TestBonds(unittest.TestCase):
  def testAddDelete(self):
    a0=molfile.Atom(name="N", resid=32)
    a1=molfile.Atom(name="C", resname="ALA", bfactor=2.5)
    
    #a0.bonds.add(a1)
    #a1.bonds.add(a0)
    a0.addbond(a1)  # also adds a0 to a1
    
    self.assertTrue( a1 in a0.bonds )
    self.assertTrue( a0 in a1.bonds )
    
    a1.delbond(a0)
    self.assertTrue( a1 not in a0.bonds )
    self.assertTrue( a0 not in a1.bonds )
    
    del a0
    del a1
    
    import gc
    gc.collect()

  def testReadWrite(self):
    a0=molfile.Atom(name="N", resid=32, anum=7)
    a1=molfile.Atom(name="C", resname="ALA", bfactor=2.5, anum=6)
    a0.addbond(a1, order=3)
    self.assertTrue( a0.getbondorder(a1)==3)
    self.assertTrue( a1.getbondorder(a0)==3)

    a1.setbondorder(a0,4)
    self.assertTrue( a0.getbondorder(a1)==4)
    self.assertTrue( a1.getbondorder(a0)==4)

    atoms=[a0,a1]

    frame=molfile.Frame(len(atoms))
    molfile.mae.write('test.mae', atoms=atoms).frame(frame).close()
    r=molfile.mae.read('test.mae')
    newatoms=r.atoms
    a0, a1 = newatoms
    self.assertTrue( a0 in a1.bonds )
    self.assertTrue( a1 in a0.bonds )
    orders=[(4,), ()]
    self.assertTrue( (r.bondorders==orders).all() )

  def tearDown(self):
    try: os.unlink('test.mae')
    except: pass

class TestPseudo(unittest.TestCase):
    def test1(self):
        mol = msys.CreateSystem()
        mol.addAtom().atomic_number=0
        mol.addAtom().atomic_number=8
        tmp = tempfile.NamedTemporaryFile(suffix='.dms')
        msys.Save(mol, tmp.name)
        r=molfile.dms.read(tmp.name)
        self.assertEqual(mol.natoms, r.natoms)

class TestFrameIterator(unittest.TestCase):

    def testMae(self):
        r=molfile.mae.read('tests/files/1vcc.mae')
        self.assertEqual(r.nframes, 1)
        self.assertEqual(len([f for f in r.frames()]), 1)
        self.assertEqual(len([f for f in r.frames()]), 1)

    def testDms(self): 
        r=molfile.dms.read('tests/files/ch4.dms')
        self.assertEqual(r.nframes, 1)
        self.assertEqual(len([f for f in r.frames()]), 1)
        self.assertEqual(len([f for f in r.frames()]), 1)
        a=r.atoms
        a[0].chain = 'x'
        a[3].chain = 'y'
        a[4].chain = 'z'
        a[0].segid = 'J'
        a[3].segid = 'G'
        dst = tempfile.NamedTemporaryFile(suffix='.dms')
        molfile.dms.write(dst.name, a).frame(r.frame(0)).close()
        r=molfile.dms.read(dst.name)
        a=r.atoms
        self.assertEqual(a[0].chain, 'x')
        self.assertEqual(a[3].chain, 'y')
        self.assertEqual(a[4].chain, 'z')
        self.assertEqual(a[0].segid, 'J')
        self.assertEqual(a[3].segid, 'G')

    def testDcd(self): 
        r=molfile.dcd.read('tests/files/alanin.dcd')
        n=len([f for f in r.frames()])
        self.assertEqual(len([f for f in r.frames()]), n)

class TestDoublePrecision(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        tmpdir = tempfile.mkdtemp()
        cls.FDTR = os.path.join(tmpdir, 'mfdbl_f.dtr')
        cls.DDTR = os.path.join(tmpdir, 'mfdbl_d.dtr')

    def tearDown(self):
        #SH.rmtree(self.FDTR, ignore_errors=True)
        #SH.rmtree(self.DDTR, ignore_errors=True)
        pass

    def check_frame(self, f, dbl, vel):
        self.assertEqual(f.fpos is None, dbl)
        self.assertEqual(f.dpos is None, not dbl)
        self.assertEqual(f.fvel is None, dbl or not vel)
        self.assertEqual(f.dvel is None, not dbl or not vel)

    def testDcd(self):
        with self.assertRaises(RuntimeError):
            r=molfile.dcd.read('tests/files/alanin.dcd', double_precision=True) 

        r=molfile.dcd.read('tests/files/alanin.dcd', double_precision=False) 
        f=next(r.frames())
        self.assertFalse(f.pos is None)
        self.assertTrue(f.vel is None)
        self.assertTrue(f.dpos is None)
        self.assertTrue(f.dvel is None)

    def testMaeDms(self):
        for dbl in True, False:
            for path, plugin in (
                    ('tests/files/ch4.dms',   molfile.dms),
                    ('tests/files/1vcc.mae',  molfile.mae),
                    ):
                r=plugin.read(path, double_precision=dbl)
                vel=r.has_velocities
                #print "path %s vel %s" % (path, vel)
                f=r.frame(0)
                self.check_frame(f, dbl, vel)
                #print path, f.fpos, f.dpos

    def testDtr(self):
        fpos = numpy.arange(3, dtype='f')
        dpos = numpy.arange(3, dtype='d')
        FDTR = molfile.dtr.write(self.FDTR, natoms=1)
        DDTR = molfile.dtr.write(self.DDTR, natoms=1)
        ff = molfile.Frame(natoms=1, double_precision=False)
        ff.pos[0][0]=0.1
        df = molfile.Frame(natoms=1, double_precision=True)
        df.dpos[0][0]=0.1
        ff.time=0; FDTR.frame(ff)
        ff.time=1; FDTR.frame(ff)
        df.time=0; DDTR.frame(df)
        df.time=1; DDTR.frame(df)
        FDTR.close()
        DDTR.close()

        for dbl in True, False:
            for path in self.FDTR, self.DDTR:
                r=molfile.dtr.read(path, double_precision=dbl)
                vel = r.has_velocities
                f = r.frame(0)
                self.check_frame(f, dbl, vel)
                # ensure no loss of precision
                if dbl and path==self.DDTR:
                    self.assertEqual(f.dpos[0][0], 0.1)
                # run a control
                if dbl and path==self.FDTR:
                    self.assertNotEqual(f.dpos[0][0], 0.1)
                    self.assertEqual(f.dpos[0][0], numpy.float32(0.1))

class TestMae(unittest.TestCase):

  def tearDown(self): 
    if os.path.isfile('out.mae'): 
      os.unlink('out.mae')
    if os.path.isfile('out.dms'): 
      os.unlink('out.dms')
 
  def testNoBonds(self):
    R=molfile.pdb.read('tests/files/h2o.pdb')
    with tempfile.NamedTemporaryFile(suffix='.mae') as tmp:
        molfile.mae.write(tmp.name, atoms=R.atoms).frame(next(R.frames())).close()
        molfile.mae.read(tmp.name)

  def testSinglePrecision(self):
      ''' mae files should preserve single precision '''
      # note - no double precision test until we have a way of writing
      # with double precision, which we don't have.
      r=molfile.dms.read('tests/files/ch4.dms')
      f=r.frame(0)
      pos=f.pos
      w=molfile.mae.write('out.mae', atoms=r.atoms)
      w.frame(r.frame(0))
      w.close()
      r2=molfile.mae.read('out.mae')
      f2=r2.frame(0)
      pos2=f2.pos
      self.assertTrue((pos==pos2).all())

  def testWithBonds(self):
    R=molfile.mae.read('tests/files/1vcc.mae')
    from time import time
    s=time()
    top=R.topology
    t=time()
    #print "topology: %f seconds" % (t-s)
    #print top
    s=time()
    order=R.bondorders
    t=time()
    #print "bondorder: %f seconds" % (t-s)
    #print order
    with tempfile.NamedTemporaryFile(suffix='.mae') as tmp:
        molfile.mae.write(tmp.name, atoms=R.atoms).frame(next(R.frames())).close()
        molfile.mae.read(tmp.name)

  def testVirtuals(self):
    R=molfile.mae.read('tests/files/tip5p.mae')
    with tempfile.NamedTemporaryFile(suffix='.mae') as tmp:
        molfile.mae.write(tmp.name, atoms=R.atoms).frame(next(R.frames())).close()
        molfile.mae.read(tmp.name)
    
  def testMultipleCt(self):
    R=molfile.mae.read('tests/files/small.mae')
    with tempfile.NamedTemporaryFile(suffix='.mae') as tmp:
        molfile.mae.write(tmp.name, atoms=R.atoms).frame(next(R.frames())).close()
        molfile.mae.read(tmp.name)

  def testPseudos(self):
    atoms=[molfile.Atom(resid=i,anum=6) for i in range(8)]
    for i in range(4): atoms[i].insertion='A'
    atoms[2].anum=0
    atoms[5].anum=0
    atoms[6].anum=0
    
    atoms[0].addbond(atoms[1])
    atoms[0].addbond(atoms[3])
    atoms[4].addbond(atoms[6])
    atoms[4].addbond(atoms[7])
    
    frame=molfile.Frame(8)
    molfile.mae.write('out.mae', atoms=atoms).frame(frame).close()
    
    # mae files reorder the particles to put pseudos in contiguous blocks
    newatoms=molfile.mae.read('out.mae').atoms
    
    old=dict([(a.resid, a.anum) for a in atoms])
    new=dict([(a.resid, a.anum) for a in newatoms])
    self.assertTrue( old==new )

class TestFrame(unittest.TestCase):
  ''' read a dcd, write out a dtr, read the dtr back in and check that
  the values are the same as the original.
  '''
  def setUp(self):
    dcd=molfile.dcd.read('tests/files/alanin.dcd')
    self.dtrpath = tempfile.mkdtemp(suffix='.dtr')
    dtr=molfile.dtr.write(self.dtrpath, natoms=dcd.natoms)
    self.x=[]
    for i, frame in enumerate(dcd.frames()):
      frame.time=i
      self.x.append(frame.pos)
      dtr.frame(frame)
    dtr.close()

  def tearDown(self):
    import shutil
    shutil.rmtree(self.dtrpath)

  def testScalars(self):
      dtr=molfile.dtr.read(self.dtrpath)
      attrs=( 'total_energy', 
              'potential_energy', 
              'kinetic_energy', 
              'extended_energy', 
              'temperature', 
              'pressure',
              )
      f=dtr.frame(0)
      for i, a in enumerate(attrs):
          self.assertEqual(getattr(f, a), None)
          setattr(f, a, i+1)
      molfile.dtr.write(self.dtrpath, natoms=dtr.natoms).frame(f).close()
      dtr=molfile.dtr.read(self.dtrpath)
      f=dtr.frame(0)
      for i, a in enumerate(attrs):
          self.assertEqual(getattr(f, a), i+1)

  def testFrame(self):
    dtr=molfile.dtr.read(self.dtrpath)
    self.assertTrue( dtr.nframes == len(self.x) )
    for i in range(dtr.nframes):
      self.assertTrue( (dtr.frame(i).pos==self.x[i]).all() )

  def testFindFrame(self):
    dtr=molfile.dtr.read(self.dtrpath)
    times=dtr.times
    self.assertTrue((times==numpy.arange(dtr.nframes)).all())

    for t,near,gt,ge,lt,le in zip(
            (32, 32.5, 0, 99, 100, -10), # time
            (32, 32,   0, 99, 99,    0), # near
            (33, 33,   1, -1, -1,    0), # gt
            (32, 33,   0, 99, -1,    0), # ge
            (31, 32,  -1, 98, 99,   -1), # lt,
            (32, 32,   0, 99, 99,   -1), # le
            ):
        self.assertEqual(findframe.at_time_near(times, t), near)
        self.assertEqual(findframe.at_time_gt(times, t), gt)
        self.assertEqual(findframe.at_time_ge(times, t), ge)
        self.assertEqual(findframe.at_time_lt(times, t), lt)
        self.assertEqual(findframe.at_time_le(times, t), le)

  def testReaderAtTime(self):
    dtr=molfile.dtr.read(self.dtrpath)

    def checkFrame(frame, time):
        if time is None:
            self.assertEqual(frame, None)
        else:
            self.assertEqual(frame.time, time)

    for t,near,gt,ge,lt,le in zip(
            (32, 32.5, 0, 99,   100, -10), # time
            (32, 32,   0, 99,   99,    0), # near
            (33, 33,   1, None, None,  0), # gt
            (32, 33,   0, 99,   None,  0), # ge
            (31, 32,None, 98,   99, None), # lt,
            (32, 32,   0, 99,   99, None), # le
            ):
        checkFrame(dtr.at_time_near(t), near)
        checkFrame(dtr.at_time_gt(t), gt)
        checkFrame(dtr.at_time_ge(t), ge)
        checkFrame(dtr.at_time_lt(t), lt)
        checkFrame(dtr.at_time_le(t), le)

class TestFrameCornerCases(unittest.TestCase):
    
    def testReaderWithNoFrameSupport(self):
        pdb=molfile.pdb.read('tests/files/h2o.pdb')
        with self.assertRaises(RuntimeError):
            pdb.at_time_near(0)

        with self.assertRaises(RuntimeError):
            pdb.times

    def testReaderWithNoTimeSupport(self):
        ''' I currently have no plugins with random access to frames but
        no implementation of read_times. '''
        pass

        pdb=molfile.pdb.read('tests/files/h2o.pdb')
        with self.assertRaises(TypeError):
            pdb.at_time_near('foobar')

    def testKeyvals(self):
        dtr=molfile.DtrReader('tests/files/force_groups.dtr')
        for i in range(dtr.nframes):
            kv = dtr.keyvals(i)
            self.assertEqual(len(list(kv.keys())), 45)
            want=[170, 169, 2678, 170, 169, 2687, 173, 169, 2678, 173, 169, 2687, 168, 169, 2687, 169, 2687, 2688, 169, 2687, 2689, 1110, 1111, 2679, 2688, 2687, 2689, 168, 169, 2678, 1110, 1111, 1112]
            self.assertEqual(
                    kv['gid__alchemical_angle_harm'].tolist(),
                    want)

    def testKeyvals2(self):
        dtr=molfile.DtrReader('tests/files/ch4.dtr')
        for i in range(dtr.nframes):
            k1 = dtr.keyvals(i)
            k2 = dict()
            dtr.frame(i, keyvals=k2)
            self.assertEqual(list(k1.keys()),list(k2.keys()))
            for k in k1:
                v1 = k1[k]
                v2 = k2[k]
                if not isinstance(v1,str):
                    v1 = v1.tolist()
                    v2 = v2.tolist()
                self.assertEqual(v1,v2)


class TestFrame2(unittest.TestCase):
  def testFrames(self):
    h=molfile.pdb.read('tests/files/h2o.pdb')
    self.assertTrue( h.natoms == 2304 )
    atoms=h.atoms
    
    for index, ts in enumerate(h.frames()):
      # FIXME: try this with non-machine representable values
      ts.pos[1] = [1,2,3]
      self.assertTrue( (ts.pos[1] == [1,2,3]).all() )
    
    for index, ts in enumerate(molfile.dcd.read('tests/files/alanin.dcd').frames()):
      pos=ts.pos
      t=index * 0.5
      ts.time = t
      self.assertTrue( ts.time == t )
    
    # test lifetime issues
    del ts
    del pos

  def testWriting(self):
    
    ### writing interface
    input=molfile.mae.read('tests/files/1vcc.mae')
    atoms=input.atoms
    molfile.psf.write('output.psf', atoms=atoms[:10]).close()
    
    
    vcc=molfile.mae.read('tests/files/1vcc.mae')
    atoms=vcc.atoms
    molfile.psf.write('1vcc.psf', atoms[:5]).close()
    
    atoms=vcc.atoms
    
    f=next(vcc.frames())
    f.box[:]=numpy.eye(3)*57
    # FIXME: need some assertions here!
    molfile.pdb.write('out.pdb', atoms[:1]).frame(f.select([3])).close()
    os.unlink('out.pdb')
    # FIMXE: here too!
    
    molfile.mae.write('my1vcc.mae', atoms=atoms).frame(f).close()
    
    inds=[1,2,3,10,20,30,40,50]
    f5=f.select(inds)
    self.assertTrue( (f.pos[inds] == f5.pos).all() )
    
    # create an atom and some frames directly
    natoms=10
    atoms=[molfile.Atom(resid=x, anum=6) for x in range(natoms)]
    frame=molfile.Frame(natoms)
    atoms[3].anum = 0
    frame.pos[:,0] = numpy.arange(natoms)
    try:
      molfile.mae.write('10.mae', atoms).frame(frame).close()
    finally:
      os.unlink('10.mae')

    os.unlink('1vcc.psf')
    os.unlink('output.psf')
    
    my=molfile.mae.read('my1vcc.mae')
    os.unlink('my1vcc.mae')
    f=next(my.frames())
    molfile.mae.write('out.mae', atoms=my.atoms).frame(f).close()

    # FIXME: check something!

class TestDtrWriter(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.PATH = os.path.join(tempfile.mkdtemp(), 'molfile_test_writer.dtr')

    def tearDown(self):
        SH.rmtree(self.PATH, ignore_errors=True)

    def testClose(self):
        d = molfile.dtr.write(self.PATH, natoms=10)
        d.close()
        del d
        d = None

    def testMixedAtomCount(self):
        d = molfile.dtr.write(self.PATH, natoms=10)
        f = molfile.Frame(20)
        f.time = 0
        with self.assertRaises(RuntimeError):
            d.frame(f)

        f = molfile.Frame(10)
        f.time = 0
        d.frame(f)

        f = molfile.Frame(20)
        f.time = 1
        with self.assertRaises(RuntimeError):
            d.frame(f)

    def testClobber(self):
        frame=molfile.Frame(10)
        molfile.dtr.write(self.PATH, natoms=10).frame(frame)
        molfile.dtr.write(self.PATH, natoms=10).frame(frame)
        molfile.dtr_clobber.write(self.PATH, natoms=10).frame(frame)
        molfile.dtr_clobber.write(self.PATH, natoms=10).frame(frame)
        molfile.dtr.write(self.PATH, natoms=10).frame(frame)

    def testCanSync(self):
        self.assertTrue(molfile.dtr.write(self.PATH, natoms=10).can_sync)

    def testAppend(self):
        # write at least 2x frames per file
        frame=molfile.Frame(10)
        dtr=molfile.dtr_append.write(self.PATH, natoms=10)
        for i in range(34):
            frame.time = i
            frame.pos[2] = i
            dtr.frame(frame)
        dtr.sync()
        dtr = molfile.dtr.read(self.PATH)
        self.assertEqual(molfile.dtr.read(self.PATH).nframes, i+1)

    def testTruncate(self):
        self.testAppend()
        dtr = molfile.dtr_append.write(self.PATH, natoms=10)
        dtr.truncate(17)
        dtr = molfile.dtr.read(self.PATH)
        self.assertEqual(molfile.dtr.read(self.PATH).nframes, 18)


    def testNoClobber(self):
        frame=molfile.Frame(10)
        self.assertFalse(os.path.exists(self.PATH))
        molfile.dtr_noclobber.write(self.PATH, natoms=10).frame(frame)
        with self.assertRaises(RuntimeError):
            molfile.dtr_noclobber.write(self.PATH, natoms=10).frame(frame)

    def testMetadata(self):
        path="/f/a/jobstep/18771855/0/diagnose.atr"
        reader=molfile.DtrReader(path)
        writer=molfile.DtrWriter(self.PATH, natoms=reader.natoms, mode=0, frames_per_file=100, metadata=reader.metadata)
        for i in range(reader.nframes):
            keyvals = reader.keyvals(i)
            writer.append(float(keyvals["CHEMICALTIME"][0]), keyvals)
        writer.sync()
        writer.close()
        r = msys.molfile.DtrReader(self.PATH)
        # vel is None when reciprocal masses are not written to the metadata file
        self.assertFalse(r.frame(0).vel is None)

class TestDtrWriterEtr(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.PATH = os.path.join(tempfile.mkdtemp(), 'molfile_test_writer.etr')

    def tearDown(self):
        SH.rmtree(self.PATH, ignore_errors=True)

    def testEtrFormat(self):
        keyvals = {'A':'B', 'X':numpy.ones(10), 'Y':numpy.array([1.0])}
        writer = msys.molfile.DtrWriter(self.PATH, 0, format=msys.molfile.DtrWriter.ETR)
        writer.append(1.0, keyvals)
        writer.append(2.0, keyvals)
        writer.append(3.0, keyvals)
        writer.sync()
        reader = msys.molfile.DtrReader(self.PATH)
        kv1 = {}
        reader.frame(0, keyvals=kv1)
        kv2 = reader.keyvals(0)
        assert(numpy.all(kv1[k] == kv2[k]) for k in kv1)
        assert(set(kv1.keys()) == set(kv2.keys()))
        assert(kv1['FORMAT'] == 'ETR_V1')
        assert(kv1['Y'][0] == 1.0)
        assert sorted(reader.metadata.keys()) == sorted(['FORMAT'] + list(keyvals.keys()))


class TestQuantizedTime(unittest.TestCase):
    def test_6659382(self):
        # they had a golden opportunity to call it ene.stk... :P
        ene = molfile.DtrReader('/f/a/job/6659382/0/energy.stk').times()
        run = molfile.DtrReader('/f/a/job/6659382/0/run.stk').times()
        self.assertEqual(ene[0], 6.48)
        self.assertEqual(run[0], 324)
        self.assertTrue(set(ene).issuperset(run))

    # too slow...   ..................................  \/
    @unittest.skipIf(os.getenv('DESRES_LOCATION')!='EN_XXX', 'Runs only from EN location')
    def test_uneven_intervals(self):
        stk = molfile.DtrReader('tests/files/uneven_intervals.stk')
        times = stk.times()
        deltas = numpy.diff(times)
        small = deltas[deltas < 0.1]
        self.assertTrue(len(small) == 404)
        self.assertTrue(deltas.min() > 0.007499)

class TestDtrFrame(unittest.TestCase):
    def test_frame_bytes_round_trip(self):
        dtr = 'tests/files/ch4.dtr'
        with open('%s/frame000000000' % dtr, 'rb') as fp:
            data = fp.read()
        d = molfile.dtr_frame_from_bytes(data)
        pos = molfile.DtrReader(dtr).frame(0).pos.flatten()
        self.assertEqual(d['POSITION'].tolist(), pos.tolist())

        for use_padding in True, False:
            data2 = molfile.dtr_frame_as_bytes(d, use_padding=use_padding)
            d2 = molfile.dtr_frame_from_bytes(data2)
            self.assertEqual(d.keys(), d2.keys())
            for k,v in d.items():
                v2 = d2[k]
                if isinstance(v, numpy.ndarray):
                    self.assertTrue((v==v2).all(), k)
                else:
                    self.assertTrue(v==v2, k)

    def test_compressed_pos_round_trip(self):
        pos = msys.Load('tests/files/2f4k.dms').positions.flatten().astype('f')
        kv = dict(POSITION=pos)
        data0 = molfile.dtr_frame_as_bytes(kv, precision=0)
        data2 = molfile.dtr_frame_as_bytes(kv, precision=1e-2)
        data3 = molfile.dtr_frame_as_bytes(kv, precision=1e-3)

        pos0 = molfile.dtr_frame_from_bytes(data0)['POSITION']
        pos2 = molfile.dtr_frame_from_bytes(data2)['POSITION']
        pos3 = molfile.dtr_frame_from_bytes(data3)['POSITION']

        self.assertTrue((pos==pos0).all())
        self.assertTrue(numpy.allclose(pos, pos2, atol=1e-2))
        self.assertFalse(numpy.allclose(pos, pos2, atol=1e-3))
        self.assertTrue(numpy.allclose(pos, pos3, atol=1e-3))
        self.assertFalse(numpy.allclose(pos, pos3, atol=1e-4))


if __name__=="__main__":
  unittest.main(verbosity=2)
