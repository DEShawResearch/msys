#!/usr/bin/env desres-exec
#{
# source `dirname $0`/../MODULES
# exec desres-cleanenv -e TMPDIR -m $PYTHON/bin -- python $0 "$@"
#}

import os, sys
TMPDIR=os.getenv('TMPDIR', 'objs/%s/x86_64' % os.getenv('DESRES_OS'))
sys.path.insert(0,os.path.join(TMPDIR, 'lib', 'python'))
#os.environ['DTRPLUGIN_VERBOSE']='1'

import sqlite3

from msys import molfile
from msys.molfile import findframe
import unittest
import numpy
import shutil as SH
import subprocess

'''
FIXME: This is sort of a mess.  There are tests for commonly used plugins,
the generic python interface, and the caching infrastructure all jumbled
together here.  
'''

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
        fids=range(100)
        random.shuffle(fids)
        frames=[f for f in self.r.frames()]
        for fid in fids:
            self.assertTrue((self.r.frame(fid).pos == frames[fid].pos).all())

class ZeromqTestCase(unittest.TestCase):
    @unittest.skipIf(True, 'Disabled as we we removed the zeromq plugin')
    def testIPC(self):
        path = 'ipc:///tmp/molfile_zeromq_test'
        natoms = 99
        
        def reader():
            r = molfile.zeromq.read(path)
            self.assertEqual(r.natoms, natoms)
            self.assertEqual(r.nframes, -1)
            i=0
            while True:
                f = r.next()
                if not f: break
                self.assertEqual(f.time, float(i))
                self.assertEqual(f.pos[37][1], float(i+1))
                self.assertEqual(f.box[1][2], float(i+2))
                if i==2:
                    self.assertEqual(f.vel[5][0], float(i+3))
                else:
                    self.assertEqual(f.vel[5][0], float(0))
                self.assertEqual(f.total_energy, i+4)
                self.assertEqual(f.potential_energy, i+5)
                self.assertEqual(f.kinetic_energy, i+6)
                self.assertEqual(f.extended_energy, i+7)
                self.assertEqual(f.pressure, i+8)
                self.assertEqual(f.temperature, i+9)
                self.assertEqual(f.pressure_tensor[2][1], i+10)
                self.assertEqual(f.virial_tensor[1][0], i+11)
                i += 1

        def writer():
            w = molfile.zeromq.write(path, natoms=natoms)
            for i in range(4):
                f = molfile.Frame(natoms, i==2)
                f.pos[37][1] = i+1
                f.time = i
                f.box[1][2] = i+2
                if i==2:
                    f.vel[5][0] = i+3
                f.total_energy = i+4
                f.potential_energy = i+5
                f.kinetic_energy = i+6
                f.extended_energy = i+7
                f.pressure = i+8
                f.temperature = i+9
                f.pressure_tensor[2][1] = i+10
                f.virial_tensor[1][0] = i+11
                w.frame(f)
            w.close()

        import multiprocessing as mp
        w=mp.Process(target=writer)
        w.start()
        r=mp.Process(target=reader)
        r.start()
        r.join()
        w.join()

class GuessFiletypeTestCase(unittest.TestCase):
    def testDtr(self):
        p=molfile.guess_filetype("foo.dtr")
        self.assertEqual(p.prettyname, "DESRES Trajectory, clobber")

class DtrGidTestCase(unittest.TestCase):
    PATH='/tmp/gids.dtr'

    def testSimple(self):
        natoms=32
        w=molfile.dtr.write(self.PATH, natoms=natoms)
        f=molfile.Frame(natoms, with_gids=True)
        gids=numpy.arange(natoms, dtype='int32')
        f.gid[:] = gids
        w.frame(f)
        f.gid[-1::-1]=gids
        w.frame(f)

        f=molfile.Frame(natoms, with_gids=False)
        w.frame(f)

        w.close()
        r=molfile.dtr.read(self.PATH, with_gids=True)
        self.assertEqual(r.natoms, natoms)
        self.assertEqual(r.nframes, 3)
        f=r.frame(0)
        self.assertTrue((f.gid==gids).all())
        f=r.frame(1)
        self.assertTrue((f.gid==gids[-1::-1]).all())
        f=r.frame(2)
        self.assertTrue((f.gid==[2**31-1]*natoms).all())

        r=molfile.DtrReader(self.PATH)
        self.assertEqual(r.natoms, natoms)
        self.assertEqual(r.nframes, 3)
        f=r.frame(0, with_gids=True)
        self.assertTrue((f.gid==gids).all())
        f=r.frame(1, with_gids=True)
        self.assertTrue((f.gid==gids[-1::-1]).all())
        f=r.frame(2, with_gids=True)
        self.assertTrue((f.gid==[2**31-1]*natoms).all())

class DtrTestCase(unittest.TestCase):
    PATH='/tmp/bad.dtr'

    def setUp(self):
        ''' creates a dtr with metadata, a timekeys file with just a header,
        and not_hashed/.ddparams '''
        SH.rmtree(self.PATH, ignore_errors=True)
        # get the dtr writer to use 1 frame per file
        self.natoms=10000000
        self.writer=molfile.dtr.write(self.PATH, natoms=self.natoms)

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
        file('%s/timekeys' % self.PATH, 'w').close()
        with self.assertRaises(RuntimeError): molfile.DtrReader(self.PATH)

    def testBadTimekeys(self):
        with file('%s/timekeys' % self.PATH, 'w') as fd: fd.write('hello')
        with self.assertRaises(RuntimeError): molfile.DtrReader(self.PATH)

    def testNoMeta(self):
        ''' test missing metadata frame.  DESRESCode#1525 '''
        # Need at least one frame in the dtr or metadata gets ignored 
        self.addFrame(0.0)
        os.unlink('%s/metadata' % self.PATH)
        molfile.DtrReader(self.PATH)

    def testOneFramePerFile(self):
        ''' test reading dtr with one frame per file.  DESRESCode#1525 '''
        self.addFrame(0.0)
        molfile.DtrReader(self.PATH).frame(0)
        self.addFrame(1.0)
        molfile.DtrReader(self.PATH).frame(0)

    def testTruncatedFrame(self):
        self.addFrame(0.0)
        self.addFrame(1.0)
        self.addFrame(2.0)
        r=molfile.DtrReader(self.PATH)
        file(r.fileinfo(1)[0], 'w').close()
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
        with file(path, 'w') as fd: print >> fd, '\0' * sz
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

    def testRelativePath(self):
        molfile.dtr.read('tests/files/relative.stk').frame(0)

    def testReplaceWithDifferentAtoms(self):
        ''' create an stk, read from it to create a cache entry, then 
        overwrite it with a new one that has a different number of atoms.
        There should be no issues. '''

        #os.environ['DTRPLUGIN_VERBOSE']='1'
        m=molfile.dtr.write('1.dtr', natoms=10)
        m.frame(molfile.Frame(10, False))
        m.close()
        with file('junk.stk', 'w') as f:
            print >> f, '1.dtr'

        #print "----- first read of junk.stk ------"
        r=molfile.dtr.read('junk.stk')
        self.assertEqual(r.nframes, 1)
        self.assertEqual(r.natoms, 10)
        self.assertEqual(r.has_velocities, False)
        self.assertEqual(r.frame(0).vel, None)

        #print "----- second read of junk.stk ------"
        r=molfile.dtr.read('junk.stk')
        self.assertEqual(r.nframes, 1)
        self.assertEqual(r.natoms, 10)
        self.assertEqual(r.has_velocities, False)
        self.assertEqual(r.frame(0).vel, None)

        #print "----- modify junk.stk ------"

        m=molfile.dtr.write('2.dtr', natoms=20)
        m.frame(molfile.Frame(20, True))
        m.frame(molfile.Frame(20, True))
        m.close()
        with file('junk.stk', 'w') as f:
            print >> f, '2.dtr'

        #print "----- read modified junk.stk ------"
        r=molfile.dtr.read('junk.stk')
        self.assertEqual(r.nframes, 2)
        self.assertEqual(r.natoms, 20)
        self.assertEqual(r.has_velocities, True)
        self.assertNotEqual(r.frame(0).vel, None)

        ## remove the velocities again
        m=molfile.dtr.write('3.dtr', natoms=20)
        m.frame(molfile.Frame(20, False))
        m.frame(molfile.Frame(20, False))
        m.close()
        with file('junk.stk', 'w') as f:
            print >> f, '3.dtr'
        r=molfile.dtr.read('junk.stk')
        self.assertEqual(r.nframes, 2)
        self.assertEqual(r.natoms, 20)
        self.assertEqual(r.has_velocities, False)
        self.assertEqual(r.frame(0).vel, None)


    def tearDown(self):
        for dtr in '1.dtr', '2.dtr', '3.dtr': 
            SH.rmtree(dtr, ignore_errors=True)
        try:
            os.unlink('junk.stk')
        except OSError:
            pass

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
        for i in xrange(dtr.nframes):
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
          allbytes = file(path).read()
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
        
        
        ## test dump/load functionality
        bytes=dtr.dump()
        oldtimes=dtr.times()
        self.assertTrue(len(bytes) < 5*1024)
        newdtr=molfile.DtrReader(bytes=bytes)
        newtimes=newdtr.times()
        self.assertEqual(newdtr.natoms, dtr.natoms)

        if (oldtimes!=newtimes).any():
          print "times differ on reload!"
          numtoshow=10
          for i, (old, new) in enumerate(zip(oldtimes, newtimes)):
              if old!=new:
                  print i, old, new
                  numtoshow -= 1
                  if not numtoshow: break
          assert False

        for i in frames:
            self.assertEqual(dtr.fileinfo(i), newdtr.fileinfo(i))
  
class TestBonds(unittest.TestCase):
  def testAddDelete(self):
    a0=molfile.Atom(name="N", resid=32)
    a1=molfile.Atom(name="C", resname="ALA", bfactor=2.5)
    
    #a0.bonds.add(a1)
    #a1.bonds.add(a0)
    a0.addbond(a1)  # also adds a0 to a1
    
    self.assert_( a1 in a0.bonds )
    self.assert_( a0 in a1.bonds )
    
    a1.delbond(a0)
    self.assert_( a1 not in a0.bonds )
    self.assert_( a0 not in a1.bonds )
    
    del a0
    del a1
    
    import gc
    gc.collect()

  def testReadWrite(self):
    a0=molfile.Atom(name="N", resid=32, anum=7)
    a1=molfile.Atom(name="C", resname="ALA", bfactor=2.5, anum=6)
    a0.addbond(a1, order=3)
    self.assert_( a0.getbondorder(a1)==3)
    self.assert_( a1.getbondorder(a0)==3)

    a1.setbondorder(a0,4)
    self.assert_( a0.getbondorder(a1)==4)
    self.assert_( a1.getbondorder(a0)==4)

    atoms=[a0,a1]

    frame=molfile.Frame(len(atoms))
    molfile.mae.write('test.mae', atoms=atoms).frame(frame).close()
    r=molfile.mae.read('test.mae')
    newatoms=r.atoms
    a0, a1 = newatoms
    self.assert_( a0 in a1.bonds )
    self.assert_( a1 in a0.bonds )
    orders=[(4,), ()]
    self.assert_( (r.bondorders==orders).all() )

  def tearDown(self):
    try: os.unlink('test.mae')
    except: pass


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
        molfile.dms.write('/tmp/_test.dms',a).frame(r.frame(0)).close()
        r=molfile.dms.read('/tmp/_test.dms')
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

    FDTR = '/tmp/mfdbl_f.dtr'
    DDTR = '/tmp/mfdbl_d.dtr'

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
        f=r.frames().next()
        self.assertNotEqual(f.pos, None)
        self.assertEqual(f.vel, None)
        self.assertEqual(f.dpos, None)
        self.assertEqual(f.dvel, None)

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
    if os.path.isfile('bonds.mae'): 
      os.unlink('bonds.mae')
 
  def testNoBonds(self):
    R=molfile.pdb.read('tests/files/h2o.pdb')
    molfile.mae.write('bonds.mae', atoms=R.atoms).frame(R.frames().next()).close()
    molfile.mae.read('bonds.mae')

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
    molfile.mae.write('out.mae', atoms=R.atoms).frame(R.frames().next()).close()
    molfile.mae.read('out.mae')

  def testVirtuals(self):
    R=molfile.mae.read('tests/files/tip5p.mae')
    molfile.mae.write('out.mae', atoms=R.atoms).frame(R.frames().next()).close()
    molfile.mae.read('out.mae')
    
  def testMultipleCt(self):
    R=molfile.mae.read('tests/files/small.mae')
    molfile.mae.write('out.mae', atoms=R.atoms).frame(R.frames().next()).close()
    molfile.mae.read('out.mae')

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
    self.assert_( old==new )

class TestFrame(unittest.TestCase):
  ''' read a dcd, write out a dtr, read the dtr back in and check that
  the values are the same as the original.
  '''
  def setUp(self):
    dcd=molfile.dcd.read('tests/files/alanin.dcd')
    dtr=molfile.dtr.write('foo.dtr', natoms=dcd.natoms)
    self.x=[]
    for i, frame in enumerate(dcd.frames()):
      frame.time=i
      self.x.append(frame.pos)
      dtr.frame(frame)
    dtr.close()

  def tearDown(self):
    import shutil
    shutil.rmtree('foo.dtr')

  def testScalars(self):
      dtr=molfile.dtr.read('foo.dtr')
      attrs=( 'total_energy', 
              'potential_energy', 
              'kinetic_energy', 
              'extended_energy', 
              'temperature', 
              'pressure',
              )
      f=dtr.frame(0)
      for i, a in enumerate(attrs):
          self.assertEquals(getattr(f, a), None)
          setattr(f, a, i+1)
      molfile.dtr.write('foo.dtr', natoms=dtr.natoms).frame(f).close()
      dtr=molfile.dtr.read('foo.dtr')
      f=dtr.frame(0)
      for i, a in enumerate(attrs):
          self.assertEquals(getattr(f, a), i+1)

  def testFrame(self):
    dtr=molfile.dtr.read('foo.dtr')
    self.assert_( dtr.nframes == len(self.x) )
    for i in range(dtr.nframes):
      self.assert_( (dtr.frame(i).pos==self.x[i]).all() )

  def testFindFrame(self):
    dtr=molfile.dtr.read('foo.dtr')
    times=dtr.times
    self.assert_((times==numpy.arange(dtr.nframes)).all())

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
    dtr=molfile.dtr.read('foo.dtr')

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

class TestFrame2(unittest.TestCase):
  def testFrames(self):
    h=molfile.pdb.read('tests/files/h2o.pdb')
    self.assert_( h.natoms == 2304 )
    atoms=h.atoms
    
    for index, ts in enumerate(h.frames()):
      # FIXME: try this with non-machine representable values
      ts.pos[1] = [1,2,3]
      self.assert_( (ts.pos[1] == [1,2,3]).all() )
    
    for index, ts in enumerate(molfile.dcd.read('tests/files/alanin.dcd').frames()):
      pos=ts.pos
      t=index * 0.5
      ts.time = t
      self.assert_( ts.time == t )
    
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
    
    f=vcc.frames().next()
    f.box[:]=numpy.eye(3)*57
    # FIXME: need some assertions here!
    molfile.pdb.write('out.pdb', atoms[:1]).frame(f.select([3])).close()
    os.unlink('out.pdb')
    # FIMXE: here too!
    
    molfile.mae.write('my1vcc.mae', atoms=atoms).frame(f).close()
    
    inds=[1,2,3,10,20,30,40,50]
    f5=f.select(inds)
    self.assert_( (f.pos[inds] == f5.pos).all() )
    
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
    f=my.frames().next()
    molfile.mae.write('out.mae', atoms=my.atoms).frame(f).close()

    # FIXME: check something!

class TestDtrWriter(unittest.TestCase):
    PATH='/tmp/molfile_test_writer.dtr'
    def tearDown(self):
        SH.rmtree(self.PATH, ignore_errors=True)

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



if __name__=="__main__":
  unittest.main(verbosity=2)