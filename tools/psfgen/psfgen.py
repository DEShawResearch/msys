import ctypes
import contextlib
import os
import msys

libpath = os.path.dirname(__file__) + '/../../libpsfgen.so'
lib=ctypes.CDLL(libpath)

''' topo_defs '''
class topo_defs(ctypes.Structure): pass
topo_defs_p = ctypes.POINTER(topo_defs)

lib.topo_defs_create.restype = topo_defs_p
lib.topo_defs_create.argtypes = []

lib.topo_defs_destroy.restype = None
lib.topo_defs_destroy.argypes = [topo_defs_p]


''' stringhash '''
class stringhash(ctypes.Structure): pass
stringhash_p = ctypes.POINTER(stringhash)

lib.stringhash_create.restype = stringhash_p

''' topo_mol '''
class topo_mol(ctypes.Structure): pass
topo_mol_p = ctypes.POINTER(topo_mol)

lib.topo_mol_create.restype = topo_mol_p

''' topo_mol_ident '''
class topo_mol_ident(ctypes.Structure):
    _fields_ = [
            ("segid", ctypes.c_char_p),
            ("resid", ctypes.c_char_p),
            ("aname", ctypes.c_char_p),
            ]

''' image_spec '''
class image_spec(ctypes.Structure):
    _fields_ = [
            ("na", ctypes.c_int),
            ("nb", ctypes.c_int),
            ("nc", ctypes.c_int),
            ("ax", ctypes.c_double),
            ("ay", ctypes.c_double),
            ("az", ctypes.c_double),
            ("bx", ctypes.c_double),
            ("by", ctypes.c_double),
            ("bz", ctypes.c_double),
            ("cx", ctypes.c_double),
            ("cy", ctypes.c_double),
            ("cz", ctypes.c_double),
            ]

''' FILE '''
class FILE(ctypes.Structure): pass
FILE_p = ctypes.POINTER(FILE)

lib.fopen.restype = FILE_p
lib.fclose.argtypes=[FILE_p]

''' callback '''
callback_p = ctypes.CFUNCTYPE(None, ctypes.c_void_p, ctypes.c_char_p)
def log(void, msg):
    print "Info:", msg
cb = callback_p(log)


''' wrapper functions '''
def read_charmm_topology(defs, path, all_caps):
    fd = lib.fopen(path, 'r')
    if not fd:
        raise RuntimeError, "Error opening charmm topology file at %s" % path
    rc = lib.charmm_parse_topo_defs(defs, fd, 1, None, cb)
    lib.fclose(fd)
    return rc

def read_pdb_residues(mol, path, alias, all_caps):
    fd = lib.fopen(path, 'r')
    if not fd:
        raise IOError, "Error opening pdb file for reading at %s" % path
    rc = lib.pdb_file_extract_residues(mol, fd, alias, all_caps, None, cb)
    lib.fclose(fd)
    return rc

def pdb_coordinates(mol, path, segid, alias, all_caps):
    fd = lib.fopen(path, 'r')
    if not fd:
        raise IOError, "Error opening pdb file for reading at %s" % path
    rc = lib.pdb_file_extract_coordinates(mol, fd, segid, alias, all_caps, None, cb)
    lib.fclose(fd)
    return rc

def _read_psf(mol, path):
    fd=lib.fopen(path, 'r')
    if not fd:
        raise IOError, "Error opening psf file for reading at %s" % path
    rc = lib.psf_file_extract(mol, fd, None, cb)
    lib.fclose(fd)
    return rc


def patch(mol, pres, targets):
    n=len(targets)
    idents = (topo_mol_ident * n)()
    for i,t in enumerate(targets):
        idents[i].segid, idents[i].resid = t
        idents[i].aname = None
    return lib.topo_mol_patch(mol, idents, n, pres, 0,0,0,0)

''' Public interface '''


class Psfgen(object):
    def __init__(self):
        self._alias = lib.stringhash_create()
        self._defs = lib.topo_defs_create()
        self._mol = lib.topo_mol_create(self._defs)
        self.all_caps = True

        lib.topo_defs_error_handler(self._defs, None, cb)
        lib.topo_mol_error_handler(self._mol, None, cb)

    def __del__(self):
        lib.topo_mol_destroy(self._mol)
        lib.topo_defs_destroy(self._defs)
        lib.stringhash_destroy(self._alias)

    def readTopology(self, path):
        ''' read a charmm topology file from the given path '''
        cb(None, "Reading charmm topology from %s" % path)
        assert 0==read_charmm_topology(self._defs, path, self.all_caps)
        lib.topo_defs_add_topofile(self._defs, path)

    @contextlib.contextmanager
    def newSegment(self, name):
        ''' begin constructing a segment with given name.  Usage:
        psf=Psfgen()
        with psf.segment('PRO') as s:
            s.pdb(...)
            s.residue(...)
        '''

        class Segment(object):
            def __init__(self, mol, name):
                self.mol = mol 
                self.name = name
            def __str__(self):
                return self.name
            def patchFirst(self, pres):
                pres=str(pres)
                print "setting patch for first residue to %s" % pres
                assert 0==lib.topo_mol_segment_first(self.mol._mol, pres)

            def patchLast(self, pres):
                pres=str(pres)
                print "setting patch for last residue to %s" % pres
                assert 0==lib.topo_mol_segment_last(self.mol._mol, pres)

            def addResidue(self, resid, resname, chain=""):
                assert 0==lib.topo_mol_residue(
                        self.mol._mol, str(resid), str(resname), str(chain))

            def mutateResidue(self, resid, resname):
                assert 0==lib.topo_mol_mutate(
                        self.mol._mol, str(resid), str(resname))

            def readSequence(self, pdb):
                path=str(pdb)
                print "reading residue sequence from pdb file %s" % path
                assert 0==read_pdb_residues(
                        self.mol._mol, path, self.mol._alias, self.mol.all_caps)

            def autoAngles(self, enable=True):
                print "%s angle autogeneration" % (
                        "enabling" if enable else "disabling")
                assert 0==lib.topo_mol_segment_auto_angles(
                        self.mol._mol, enable)

            def autoDihedrals(self, enable=True):
                print "%s dihedral autogeneration" % (
                        "enabling" if enable else "disabling")
                assert 0==lib.topo_mol_segment_auto_dihedrals(
                        self.mol._mol, enable)

        assert 0==lib.topo_mol_segment(self._mol, name)
        yield Segment(self, name)
        assert 0==lib.topo_mol_end(self._mol)

    def guessCoordinates(self):
        assert 0==lib.topo_mol_guess_xyz(self._mol)

    def regenerateAngles(self):
        assert 0==lib.topo_mol_regenerate_angles(self._mol)

    def regenerateDihedrals(self):
        assert 0==lib.topo_mol_regenerate_dihedrals(self._mol)

    def aliasResidue(self, **kwds):
        for old,v in kwds.items():
            new=str(v)
            print "aliasing residue %s to %s" % (old, new)
            assert 0==lib.extract_alias_residue_define(self._alias, old, new)

    def aliasAtom(self, residue, **kwds):
        residue=str(residue)
        for old,v in kwds.items():
            new=str(v)
            print "aliasing residue %s atom %s to %s" % (residue, old, new)
            assert 0==lib.extract_alias_atom_define(self._alias, residue, old, new)

    def readPDB(self, path, segid=None):
        path=str(path)
        print "reading coordinates from pdb file %s" % path
        if segid is not None:
            segid = str(segid)
        assert 0==pdb_coordinates(
                self._mol, str(path), segid, self._alias, self.all_caps)

    def readPSF(self, path):
        path=str(path)
        print "Reading structure from psf file %s" % path
        assert 0==_read_psf(self._mol, path)

    def writePSF(self, path):
        path=str(path)
        fd=lib.fopen(path, "w")
        if fd is None:
            raise IOError, "Error opening psf file for writing at %s" % path
        print "Writing psf file to %s" % path
        rc=lib.topo_mol_write_psf(self._mol, fd, 0,0, None, cb)
        lib.fclose(fd)
        assert rc==0
        print "psf file complete"

    def writePDB(self, path):
        path=str(path)
        fd=lib.fopen(path, "w")
        if fd is None:
            raise IOError, "Error opening pdb file for writing at %s" % path
        print "Writing pdb file to %s" % path
        rc=lib.topo_mol_write_pdb(self._mol, fd, None, cb)
        lib.fclose(fd)
        assert rc==0
        print "pdb file complete"

    def read(self, type, path):
        type=str(type)
        path=str(path)
        fd=lib.fopen(path, 'r')
        if not fd:
            raise IOError, "Error opening %s file for reading at %s" % (
                    type, path)
        print "Reading %s file %s" % (type,path)
        segid=None
        all_caps=0
        coordinatesonly=0
        residuesonly=0
        rc=lib.topo_mol_read_plugin(self._mol, 
                                    type, path,
                                    type, path,
                                    segid,
                                    self._alias, 
                                    all_caps,
                                    coordinatesonly, residuesonly,
                                    0, cb)
        lib.fclose(fd)
        assert 0==rc

    def write(self, type, path):
        type=str(type)
        path=str(path)
        img = image_spec()
        ctypes.memset(ctypes.addressof(img), 0, ctypes.sizeof(img))
        img.na = img.nb = img.nc = 1
        print "Writing structure using '%s' plugin to %s" % (type,path)
        assert 0==lib.topo_mol_write_plugin(self._mol, type, path,
                ctypes.byref(img), None, cb)
        print "Writing structure complete."

    def applyPatch(self, pres, *args):
        pres=str(pres)
        targets=list()
        for arg in args:
            try:
                segid, resid = map(str, arg)
            except:
                raise ValueError, "Argument '%s' is not of the form (SEGID,RESID)" % arg
            targets.append((segid,resid))
        print "applying patch %s to %d residues" % (pres, len(args))
        assert 0==patch(self._mol, pres, targets)

    def applyResidue(self, res, segid, resid):
        ''' Mutate the residue at segid, resid to the topology defined in res
        '''
        idents = (topo_mol_ident * 1)()
        idents[0].segid = str(segid)
        idents[0].resid = str(resid)
        res=str(res)
        assert 0==lib.topo_mol_patch_residue(self._mol, idents, res)

class Params(object):
    def __init__(self):
        self.bond = list()
        self.angle = list()
        self.dihedral = list()
        self.improper = list()
        self.cmap = list()
        self.nonbonded = list()

        self.convert_sigma = 2.0 / 2.0**(1./6.)

    def parse(self, fobj):
        mode=self.error
        handlers={
                'BOND' : self.parse_bond,
                'ANGL' : self.parse_angle,
                'DIHE' : self.parse_dihedral,
                'IMPR' : self.parse_improper,
                'CMAP' : self.parse_cmap,
                'NONB' : self.parse_nonbonded,
                }

        lines=fobj.readlines()
        i, n = 0, len(lines)
        while i<n:
            line = lines[i]
            i += 1
            if line.startswith('*'):
                print line,
                continue
            line=line.split('!')[0].strip()
            if not line: 
                continue
            line=line.upper()
            if line[:8]=='READ RTF':
                begin=i
                while i<n:
                    line=lines[i]
                    if line.upper().strip()=='END':
                        i += 1
                        break
                    i += 1
                end=i
                #print "skipped lines %d-%d in %s" % (begin,end,fobj.name)
                continue
            elif line[:9]=='READ PARA':
                continue
            elif line[:5]=='CUTNB':
                continue
            elif line[:5]=='HBOND':
                continue
            elif line[:3]=='END':
                continue
            start=line[:4]
            func=handlers.get(start)
            if func is not None: 
                mode=func
                continue
            try:
                mode(line)
            except:
                print "Error parsing line %d in %s:\n%s" % (i, fobj.name, line)
                raise

    def parse_bond(self, line):
        t1, t2, fc, r0 = line.split()
        params={'fc' : float(fc), 'r0' : float(r0) }
        self.bond.append(((t1,t2), params))
        self.bond.append(((t2,t1), params))

    def parse_angle(self, line):
        elems = line.split()
        t1, t2, t3 = elems[:3]
        fc=float(elems[3])
        theta0=float(elems[4])
        params={'fc' : fc, 'theta0' : theta0 }
        if len(elems)>5:
            ufc, ur0 = map(float, elems[5:])
            params.update(ufc=ufc, ur0=ur0)
        self.angle.append(((t1,t2,t3),params))
        self.angle.append(((t3,t2,t1),params))

    def parse_dihedral(self, line):
        t1, t2, t3, t4, fc, n, delta = line.split()
        types=(t1,t2,t3,t4)
        fc=float(fc)
        absfc=fc
        n=int(n)
        delta=float(delta)
        if delta==0:
            phi0=0.0
        elif delta==180:
            phi0=0.0
            fc=-fc
        else:
            phi0=delta
        keys='phi0', 'fc0', 'fc1', 'fc2', 'fc3', 'fc4', 'fc5', 'fc6'
        prev=() if not len(self.dihedral) else self.dihedral[-1]
        if prev and prev[0]==types and prev[1]['phi0']==phi0:
            params=prev[1]
            params[keys[1  ]]+=absfc
            params[keys[n+1]]=fc
        else:
            params=dict()
            params[keys[0  ]]=phi0
            params[keys[1  ]]=absfc
            params[keys[n+1]]=fc
            self.dihedral.append((types, params))

    def parse_improper(self, line):
        t1, t2, t3, t4, fc, _, phi0 = line.split()
        types=(t1,t2,t3,t4)
        fc=float(fc)
        phi0=float(phi0)
        params={'fc' : fc, 'phi0' : phi0 }
        self.improper.append((types, params))

    def parse_cmap(self, line):
        elems=line.split()
        if len(elems)==9:
            self.cmap.append((elems[:-1],[]))
        else:
            self.cmap[-1][-1].extend(map(float,elems))

    def parse_nonbonded(self, line):
        elems=line.split()
        key=elems[0]
        ignored=elems[1]
        epsilon=-float(elems[2])
        sigma=self.convert_sigma * float(elems[3])
        params={'sigma' : sigma, 'epsilon' : epsilon}
        if len(elems)>4:
            _, eps14, rmin14 = elems[4:]
            eps14=-float(eps14)
            sig14=self.convert_sigma * float(rmin14)
            params['sigma14'] = sig14
            params['epsilon14'] = eps14
        self.nonbonded.append((key, params))

    def error(self, line):
        print "error!", line

def assign(params, mol):
    nb=mol.table('nonbonded')
    types=[p['type'] for p in nb.params.params]
    atypes=[types[nb.term(i).param.id] for i in range(mol.natoms)]

    from time import time

    T0=time()
    # parameterize bonds
    stretch=mol.addTableFromSchema('stretch_harm')
    plist=[]
    for t, d in params.bond:
        p = stretch.params.addParam()
        p['fc'] = d['fc']
        p['r0'] = d['r0']
        plist.append(p)

    for b in mol.bonds:
        ai, aj = b.atoms
        t = (atypes[ai.id], atypes[aj.id])
        for i, p in enumerate(params.bond):
            if p[0] == t:
                atoms = [ai,aj] if ai.id < aj.id else [aj,ai]
                stretch.addTerm(atoms, plist[i])
                break
        else:
            raise RuntimeError, "no match for atoms %s" % (b.atoms,)

    T1=time()
    # parameterize angles and urey-bradley terms
    angles=mol.addTableFromSchema('angle_harm')
    plist=[]
    ulist=dict()
    for t, d in params.angle:
        p = angles.params.addParam()
        p['fc'] = d['fc']
        p['theta0'] = d['theta0']
        plist.append(p)
        if 'ufc' in d:
            p = stretch.params.addParam()
            p['fc']=d['ufc']
            p['r0']=d['ur0']
            ulist[len(plist)-1]=p

    terms=[]
    keyhash=dict()
    if 'angles' in mol.auxtable_names:
        terms=mol.auxtable('angles').params
    for a in terms:
        ids=[a[x] for x in 'p0', 'p1', 'p2']
        if ids[0]>ids[2]:
            ids[0], ids[2] = ids[2], ids[0]
        t=tuple(atypes[i] for i in ids)
        p = keyhash.get(t)
        if p is not None:
            i,p = p
            angles.addTerm([mol.atom(x) for x in ids], plist[i])
            if 'ufc' in p[1]:
                atoms=[mol.atom(ids[x]) for x in 0,2]
                stretch.addTerm(atoms, ulist[i])
            continue
        for i, p in enumerate(params.angle):
            if p[0] == t:
                angles.addTerm([mol.atom(x) for x in ids], plist[i])
                if 'ufc' in p[1]:
                    atoms=[mol.atom(ids[x]) for x in 0,2]
                    stretch.addTerm(atoms, ulist[i])
                keyhash[t]=(i,p)
                break
        else:
            raise RuntimeError, "no match for atoms %s with types %s" % (
                    ids, t)

    T2=time()
    # dihedrals
    trig=mol.addTableFromSchema('dihedral_trig')
    for _,d in params.dihedral:
        p = trig.params.addParam()
        for k,v in d.items(): p[k]=v
    terms=list()
    keyhash=dict()
    if 'dihedrals' in mol.auxtable_names:
        terms=mol.auxtable('dihedrals').params
    for a in terms:
        ids=[a[x] for x in 'p0', 'p1', 'p2', 'p3']
        if ids[0]>ids[3]:
            ids[0], ids[3] = ids[3], ids[0]
        atoms=[mol.atom(x) for x in ids]
        t=tuple(atypes[i] for i in ids)
        p=keyhash.get(t)
        if p is not None:
            trig.addTerm(atoms, p)
            continue
        tr=tuple(reversed(t))
        for i,(key,d) in enumerate(params.dihedral):
            wild=key[0]=='X' and key[3]=='X'
            mid=key[1:3]
            if t==key or tr==key or (wild and (t[1:3]==mid or tr[1:3]==mid)):
                p=trig.params.param(i)
                keyhash[t]=p
                keyhash[tr]=p
                trig.addTerm(atoms, p)
                break
        else:
            raise RuntimeError, "no match for dihedral %s with types %s" % (ids,t)

    T3=time()
    # impropers
    impr=mol.addTableFromSchema('improper_harm')
    for _,d in params.improper:
        p = impr.params.addParam()
        for k,v in d.items(): p[k]=v
    terms=list()
    if 'impropers' in mol.auxtable_names:
        terms=mol.auxtable('impropers').params
    for a in terms:
        ids=[a[x] for x in 'p0', 'p1', 'p2', 'p3']
        if ids[0]>ids[3]:
            ids[0], ids[3] = ids[3], ids[0]
        t=tuple(atypes[i] for i in ids)
        tr=tuple(reversed(t))
        for i,(key,d) in enumerate(params.improper):
            wild=key[1]=='X' and key[2]=='X'
            out=key[0::3]
            if t==key or tr==key or (wild and (t[0::3]==out or tr[0::3]==out)):
                atoms=[mol.atom(x) for x in ids]
                impr.addTerm(atoms, impr.params.param(i))
                break
        else:
            raise RuntimeError, "no match for improper %s with types %s" % (ids,t)

    T4=time()
    # exclusions: assuming 1-2, 1-3, and 1-4
    excl=mol.addTableFromSchema('exclusion')
    s23=set()
    s14=set()
    for table in 'stretch_harm', 'angle_harm', 'dihedral_trig':
        s=s14 if mol.table(table).natoms==4 else s23
        for t in mol.table(table).terms:
            atoms=t.atoms
            ai=atoms[0]
            aj=atoms[-1]
            s.add(tuple(sorted((ai.id, aj.id))))
    for ai,aj in s23.union(s14):
        excl.addTerm([mol.atom(ai), mol.atom(aj)], None)

    T5=time()
    # parameterize nonbonded
    nb=mol.addNonbondedFromSchema('vdw_12_6', 'arithmetic/geometric')
    nb.params.addProp('sigma', float)
    nb.params.addProp('epsilon', float)
    nbdict=dict()
    for p in nb.params.params:
        t = p['type']
        for key, d in params.nonbonded:
            if t==key:
                p['sigma'] = d['sigma']
                p['epsilon'] = d['epsilon']
                nbdict[t]=d
                break
        else:
            raise RuntimeError, "No nonbonded param for '%s'" % t

    T6=time()

    # pairs: assuming 1-4 pairs with full scaling
    pairs=mol.addTableFromSchema('pair_12_6_es')
    keyhash=dict()
    for a1, a2 in s14.difference(s23):
        t1=atypes[a1]
        t2=atypes[a2]
        sigeps = keyhash.get((t1,t2))
        if sigeps is None:
            d1=nbdict[t1]
            d2=nbdict[t2]
            s1=d1.get('sigma14', d1['sigma'])
            s2=d2.get('sigma14', d2['sigma'])
            e1=d1.get('epsilon14', d1['epsilon'])
            e2=d2.get('epsilon14', d2['epsilon'])
            sij=0.5*(s1+s2)
            eij=(e1*e2)**0.5
            sigeps = (sij,eij)
            keyhash[(t1,t2)]=sigeps
            keyhash[(t2,t1)]=sigeps
        else:
            sij, eij = sigeps
        p=pairs.params.addParam()
        ai=mol.atom(a1)
        aj=mol.atom(a2)
        p['aij'] = sij**12 * eij * 4.0
        p['bij'] = sij**6  * eij * 4.0
        p['qij'] = ai.charge * aj.charge
        pairs.addTerm([ai,aj], p)

    T7=time()

    # cmaps
    cmaps=mol.addTableFromSchema('torsiontorsion_cmap')
    for i in range(6):
        p=cmaps.params.addParam()
        p['cmapid']="cmap%d" % (i+1)

    terms=list()
    used_tables=set()
    if 'cmaps' in mol.auxtable_names:
        terms=mol.auxtable('cmaps').params
    for a in terms:
        ids=[a[x] for x in 'p0', 'p1', 'p2', 'p3', 'p4', 'p5', 'p6', 'p7']
        t=list(atypes[i] for i in ids)
        for i, (key, d) in enumerate(params.cmap):
            if t==key:
                used_tables.add(i)
                atoms=list(mol.atom(x) for x in ids)
                param=cmaps.params.param(i)
                cmaps.addTerm(atoms, param)
                break
        else:
            raise RuntimeError, "No cmap for '%s'" % t

    for i in used_tables:
        key, d = params.cmap[i]
        t=msys.CreateParamTable()
        t.addProp('phi', float)
        t.addProp('psi', float)
        t.addProp('energy', float)
        j=0
        for phi in range(-180, 180,15):
            for psi in range(-180,180,15):
                p=t.addParam()
                p['phi']=phi
                p['psi']=psi
                p['energy']=d[j]
                j += 1
        mol.addAuxTable('cmap%d' % (i+1), t)


    #print "bonds        %12f" % (T1-T0)
    #print "angles       %12f" % (T2-T1)
    #print "dihedrals    %12f" % (T3-T2)
    #print "impropers    %12f" % (T4-T3)
    #print "exclusions   %12f" % (T5-T4)
    #print "nonbonded    %12f" % (T6-T5)
    #print "pairs        %12f" % (T7-T6)

    for t in 'angles', 'dihedrals', 'impropers', 'cmaps':
        if t in mol.auxtable_names:
            mol.delAuxTable(t)

Params.assign = assign

def AssignDMS(ifile, ofile, *paramfiles):
    print "Reading DMS file from %s" % ifile
    mol=msys.LoadDMS(ifile)
    params=Params()
    for p in paramfiles:
        with file(p) as fd:
            print "Parsing parameter file at", p
            params.parse(fd)
    print "Assigning parameters..."
    params.assign(mol)
    mol.coalesceTables()
    mol=mol.clone()
    print "Writing DMS file to %s" % ofile
    msys.SaveDMS(mol, ofile)
    print "Done."

