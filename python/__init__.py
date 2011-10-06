import _msys

class Handle(object):
    __slots__ = ('_ptr', '_id')

    def __init__(self, _ptr, _id):
        self._ptr = _ptr
        self._id = _id

    def __eq__(self, x): return self._id==x._id and self._ptr==x._ptr
    def __ne__(self, x): return self._id!=x._id or  self._ptr!=x._ptr

    @property
    def id(self): return self._id
    
    @property
    def system(self): return System(self._ptr)


class Bond(Handle):
    def destroy(self):
        self._ptr.delBond(self.id)

    @property
    def bond(self): return self._ptr.bond(self.id)

    @property
    def first(self): return Atom(self._ptr, self.bond.i)

    @property
    def second(self): return Atom(self._ptr, self.bond.j)

    @property
    def order(self): return self.bond.order
    @order.setter
    def order(self, val): self.bond.order = val

    @property
    def atoms(self): return self.first, self.second



class Atom(Handle):
    __slots__ = ()

    def destroy(self):
        self._ptr.delAtom(self.id)

    def __setitem__(self, key, val):
        p=self._ptr
        p.setAtomProp(self.id, key, val)

    def __getitem__(self, key):
        p=self._ptr
        return p.getAtomProp(self.id, key)
    
    def __contains__(self, key):
        return not _msys.bad(self._ptr.atomPropIndex(key))

    def addBond(self, other):
        assert self._ptr == other._ptr
        return Bond(self._ptr, self._ptr.addBond(self.id, other.id))

    @property
    def atom(self): return self._ptr.atom(self.id)

    @property
    def residue(self): return Residue(self._ptr, self.atom.residue)
    @residue.setter
    def residue(self, res): self._ptr.setResidue(self.id, res.id)

    @property
    def bonds(self):
        return [Bond(self._ptr, i) for i in self._ptr.bondsForAtom(self.id)]

    @property
    def bonded_atoms(self):
        return [Atom(self._ptr, i) for i in self._ptr.bondedAtoms(self.id)]


    charge = property( lambda self: self.atom.charge,
                       lambda self, val: setattr(self.atom, 'charge', val))
#for attr in (
    #'gid', 'fragid', 'x', 'y', 'z', 'charge', 'vx', 'vy', 'vz',
    #'mass', 'chargeB', 'atomic_number', 'formal_charge',
    #'moiety', 'alchemical', 'name'):
    #Atom.__dict__[attr]=property( 
            #lambda self: getattr(self.atom, attr),
            #lambda self, val: setattr(self.atom, attr, val))
    
class Residue(Handle):

    def destroy(self):
        self._ptr.delResidue(self.id)

    def addAtom(self):
        return Atom(self._ptr, self._ptr.addAtom(self.id))

    @property
    def atoms(self):
        return [Atom(self._ptr, i) for i in self._ptr.atomsForResidue(self.id)]

class Chain(Handle):

    def destroy(self):
        self._ptr.delChain(self.id)

    def addResidue(self):
        return Residue(self._ptr, self._ptr.addResidue(self.id))

    @property
    def residues(self):
        return [Residue(self._ptr, i) for i in self._ptr.residuesForChain(self.id)]

class Param(object):
    __slots__=('_ptr', '_id')

    def __init__(self, _ptr, _id):
        self._ptr = _ptr
        self._id  = _id

    def __eq__(self, x): return self._id==x._id and self._ptr==x._ptr
    def __ne__(self, x): return self._id!=x._id or  self._ptr!=x._ptr

    @property
    def id(self): return self._id
    
    @property
    def table(self): return ParamTable(self._ptr)

    def __setitem__(self, key, val):
        p=self._ptr
        p.setProp(self.id, p.propIndex(key), val)

    def __getitem__(self, key):
        p=self._ptr
        return p.getProp(self.id, p.propIndex(key))

class ParamTable(object):
    __slots__=('_ptr',)

    def __init__(self, _ptr):
        ''' Construct from ParamTablePtr.
        Do not invoke directly; use CreateParamTable() instead.  '''
        self._ptr = _ptr

    def __eq__(self, x): return self._ptr==x._ptr
    def __ne__(self, x): return self._ptr!=x._ptr

    def addParam(self):
        ''' add and return a new entry in the parameter table '''
        return Param(self._ptr, self._ptr.addParam())

    def addProp(self, name, type):
        return self._ptr.addProp(name,type)

    @property
    def props(self):
        p=self._ptr
        return [p.propName(i) for i in range(p.propCount())]

    @property
    def prop_types(self):
        p=self._ptr
        return [p.propType(i) for i in range(p.propCount())]

    def propType(self, name):
        p=self._ptr
        return p.propType(p.propIndex(name))

    def __getitem__(self, id):
        ''' fetch the Param with the given id '''
        return Param(self._ptr, id)

    def __len__(self):
        return self._ptr.paramCount()


class Term(object):
    __slots__=('_ptr', '_id')

    def __init__(self, _ptr, _id):
        self._ptr = _ptr
        self._id  = _id

    def __eq__(self, x): return self._id==x._id and self._ptr==x._ptr
    def __ne__(self, x): return self._id!=x._id or  self._ptr!=x._ptr

    @property
    def id(self): return self._id

    @property
    def param(self): 
        id=self._ptr.param(self.id)
        if _msys.bad(id): return None
        return Param(self._ptr.paramTable(), id)
    @param.setter
    def param(self, val):
        if val is None: id = None
        else: id = val.id
        self._ptr.setParam(self.id, id)

    @property
    def atoms(self):
        return [Atom(self._ptr.system(), i) for i in self._ptr.atoms(self.id)]

    @property
    def table(self): return TermTable(self._ptr)

    def __getitem__(self, attr):
        col = self._ptr.propIndex(attr)
        if _msys.bad(col):
            col = self._ptr.termPropIndex(attr)
            if _msys.bad(col):
                raise ValueError, "No such property '%s'" % attr
            return self._ptr.getTermProp(self.id, col)
        return self._ptr.getProp(self.id, col)

    def __setitem__(self, attr, val):
        col = self._ptr.propIndex(attr)
        if _msys.bad(col):
            col = self._ptr.termPropIndex(attr)
            if _msys.bad(col):
                raise ValueError, "No such property '%s'" % attr
            self._ptr.setTermProp(self.id, col, val)
        else:
            self._ptr.setProp(self.id, col, val)



class TermTable(object):
    __slots__=('_ptr',)

    def __init__(self, _ptr):
        ''' Construct from TermTablePtr.
        Do not invoke directly; use System.addTable or System.table instead '''
        self._ptr = _ptr

    def __eq__(self, x): return self._ptr==x._ptr
    def __ne__(self, x): return self._ptr!=x._ptr

    @property
    def param_table(self): return ParamTable(self._ptr.paramTable())

    @property
    def system(self): return System(self._ptr.system())

    @property
    def props(self): 
        p=self._ptr
        return [p.propName(i) for i in range(p.propCount())]

    @property
    def prop_types(self): 
        p=self._ptr
        return [p.propType(i) for i in range(p.propCount())]

    @property
    def term_props(self): 
        p=self._ptr
        return [p.termPropName(i) for i in range(p.termPropCount())]

    @property
    def term_prop_types(self): 
        p=self._ptr
        return [p.termPropType(i) for i in range(p.termPropCount())]

    @property
    def natoms(self): return self._ptr.atomCount()

    @property
    def category(self): return self._ptr.category
    @category.setter
    def category(self, val): self._ptr.category=val

    @property
    def terms(self):
        return [Term(self._ptr, i) for i in self._ptr.terms()]

    def addProp(self, name, type):
        return self._ptr.addProp(name, type)

    def propType(self, name):
        p=self._ptr
        return self._ptr.propType(p.propIndex(name))

    def termPropType(self, name):
        p=self._ptr
        return p.termPropType(p.termPropIndex(name))


    def addTermProp(self, name, type):
        return self._ptr.addTermProp(name, type)


    def addTerm(self, atoms, param = None):
        if param is not None:
            param = param.id
        n=len(atoms)
        ids=[None]*n
        for i in range(n):
            a=atoms[i]
            ids[i]=a.id
            if a.system != self.system:
                raise RuntimeError, "Cannot add atoms from different system"
        return Term(self._ptr, self._ptr.addTerm(ids, param))

class System(object):
    def __init__(self, _ptr):
        ''' Construct from SystemPtr.
        Do not invoke directly; use CreateSystem() instead.
        '''
        self._ptr = _ptr

    def __eq__(self, x): return self._ptr == x._ptr
    def __ne__(self, x): return self._ptr != x._ptr

    def addAtom(self):
        ''' add and return a new Atom in its own residue '''
        return self.addResidue().addAtom()

    def addResidue(self):
        ''' add and return a new Residue in its own chain '''
        return self.addChain().addResidue()

    def addChain(self):
        ''' add and return a new Chain '''
        return Chain(self._ptr, self._ptr.addChain())

    def delAtoms(self, elems):
        self._ptr.delAtoms(e.id for e in elems)

    def delBonds(self, elems):
        self._ptr.delBonds(e.id for e in elems)

    def delResidues(self, elems):
        self._ptr.delResidues(e.id for e in elems)

    def delChains(self, elems):
        self._ptr.delChains(e.id for e in elems)

    def atom(self, id):
        ''' return the atom with the specified id '''
        return Atom(self._ptr, id)

    @property
    def cell(self):
        return self._ptr.global_cell

    @property
    def atoms(self):
        ''' return list of all atoms in the system '''
        ptr=self._ptr
        return [Atom(ptr, i) for i in ptr.atoms()]

    @property
    def bonds(self):
        ''' return list of all bonds in the system '''
        ptr=self._ptr
        return [Bond(ptr, i) for i in ptr.bonds()]

    @property
    def residues(self):
        ''' return list of all residues in the system '''
        ptr=self._ptr
        return [Residue(ptr, i) for i in ptr.residues()]

    @property
    def chains(self):
        ''' return list of all chains in the system '''
        ptr=self._ptr
        return [Chain(ptr, i) for i in ptr.chains()]

    @property
    def residues(self):
        ''' return list of all residues in the sytem '''
        ptr=self._ptr
        return [Residue(ptr, i) for i in ptr.residues()]

    def addAtomProp(self, name, type):
        ''' add a custom atom property with the given name and type.
        type should be int, float, or str.
        '''
        return self._ptr.addAtomProp(name, type)

    @property
    def atom_props(self):
        ''' return the list of custom atom properties. '''
        p=self._ptr
        return [p.atomPropName(i) for i in range(p.atomPropCount())]

    def atomPropType(self, name):
        return self._ptr.atomPropType(self._ptr.atomPropIndex(name))

    @property
    def atom_prop_types(self):
        ''' return the types of the custom atom properties '''
        p=self._ptr
        return [p.atomPropType(i) for i in range(p.atomPropCount())]

    def addTable(self, name, natoms, params = None):
        if params is None: params = CreateParamTable()
        return TermTable(self._ptr.addTable(name, natoms, params._ptr))

    def addTableFromSchema(self, type, name = None):
        if name is None: name=type
        return TermTable(self._ptr.addTableFromSchema(type,name))

    def addNonbondedFromSchema(self, funct, rule):
        return TermTable(self._ptr.addNonbondedFromSchema(funct,rule))

    def atomselect(self, seltext):
        p=self._ptr
        ids=p.atomselect(seltext)
        return [Atom(p,i) for i in ids]

def CreateSystem():
    return System(_msys.SystemPtr.create())

def CreateParamTable():
    return ParamTable(_msys.ParamTablePtr.create())

def LoadDMS(path, structure_only = False):
    return System(_msys.ImportDMS(path, structure_only ))

def LoadMAE(path, ignore_unrecognized = False):
    return System(_msys.ImportMAE(path, ignore_unrecognized ))

def SaveDMS(system, path):
    _msys.ExportDMS(system._ptr, path)


''' customize Vec3 '''
from _msys import Vec3
def __vec3_getitem(self, key):
    tmp = [self.x, self.y, self.z]
    return tmp[key]

def __vec3_setitem(self, key, val):
    tmp = [self.x, self.y, self.z]
    tmp[key] = val
    self.x, self.y, self.z = tmp
Vec3.__getitem__ = __vec3_getitem
Vec3.__setitem__ = __vec3_setitem


