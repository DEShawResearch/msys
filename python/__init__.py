
'''
This is the high-level Python interface for msys, intended for use
by chemists.
'''

import _msys
import numpy
import sys

from _msys import GlobalCell, NonbondedInfo, version, hexversion
from _msys import RadiusForElement, MassForElement
from atomsel import Atomsel

class Handle(object):
    __slots__ = ('_ptr', '_id')

    def __init__(self, _ptr, _id):
        self._ptr = _ptr
        self._id = _id

    def __eq__(self, x): 
        return self.__class__==type(x) and self._id==x._id and self._ptr==x._ptr

    def __ne__(self, x): 
        return self.__class__!=type(x) or self._id!=x._id or self._ptr!=x._ptr

    def __hash__(self): return hash((self._ptr, self._id))

    @property
    def id(self): 
        ''' id of this object in the parent System '''
        return self._id
    
    @property
    def system(self): 
        ''' parent System '''
        return System(self._ptr)

def __add_properties(cls, *names):
    for attr in names:
        setattr(cls, attr, property(
            lambda self, x=attr: getattr(self.data(), x),
            lambda self, val, x=attr: setattr(self.data(), x, val),
            doc=attr))

class Bond(Handle):
    __slots__ = ()

    def __repr__(self): return '<Bond %d>' % self._id

    def data(self): return self._ptr.bond(self._id)

    def remove(self): 
        ''' remove this Bond from the System '''
        self._ptr.delBond(self._id)

    @property
    def first(self): 
        ''' first Atom in the bond (the one with lower id) '''
        return Atom(self._ptr, self.data().i)

    @property
    def second(self): 
        ''' second Atom in the bond (the one with higher id) '''
        return Atom(self._ptr, self.data().j)

    @property
    def atoms(self): 
        ''' Atoms in this Bond '''
        return self.first, self.second

    def __setitem__(self, key, val):
        ''' set custom Bond property '''
        self._ptr.setBondProp(self._id, key, val)

    def __getitem__(self, key):
        ''' get custom Bond property '''
        return self._ptr.getBondProp(self._id, key)
    
    def __contains__(self, key):
        ''' does custom Bond property exist? '''
        return not _msys.bad(self._ptr.bondPropIndex(key))


__add_properties(Bond, 'order')
__add_properties(Bond, 'resonant_order')

class Atom(Handle):
    __slots__ = ()

    def __repr__(self): return '<Atom %d>' % self._id

    def data(self): return self._ptr.atom(self._id)

    def remove(self):
        ''' remove this Atom from the System '''
        self._ptr.delAtom(self._id)

    def __setitem__(self, key, val):
        ''' set atom property key to val '''
        self._ptr.setAtomProp(self._id, key, val)

    def __getitem__(self, key):
        ''' get atom property key '''
        return self._ptr.getAtomProp(self._id, key)
    
    def __contains__(self, key):
        ''' does atom property key exist? '''
        return not _msys.bad(self._ptr.atomPropIndex(key))

    def addBond(self, other):
        ''' create and return a Bond from self to other '''
        assert self._ptr == other._ptr
        return Bond(self._ptr, self._ptr.addBond(self._id, other.id))

    def findBond(self, other):
        ''' Find the bond between self and Atom other; None if not found '''
        return self.system.findBond(self, other)

    @property
    def pos(self):
        ''' position '''
        return numpy.array((self.x, self.y, self.z))
    @pos.setter
    def pos(self, xyz):
        self.x, self.y, self.z = map(float, xyz)

    @property
    def vel(self):
        ''' velocity '''
        return numpy.array((self.vx, self.vy, self.vz))
    @vel.setter
    def vel(self, xyz):
        self.vx, self.vy, self.vz = map(float, xyz)

    @property
    def residue(self): return Residue(self._ptr, self.data().residue)
    @residue.setter
    def residue(self, res): self._ptr.setResidue(self._id, res.id)

    @property
    def bonds(self):
        ''' Bonds connected to this atom '''
        return [Bond(self._ptr, i) for i in self._ptr.bondsForAtom(self._id)]

    @property
    def bonded_atoms(self):
        ''' Atoms bonded to this atom '''
        return [Atom(self._ptr, i) for i in self._ptr.bondedAtoms(self._id)]

    @property
    def nbonds(self):
        ''' number of bonds to this atom '''
        return self._ptr.bondCountForAtom(self._id)


__add_properties(Atom, 
        'fragid', 'x', 'y', 'z', 'charge',
        'vx', 'vy', 'vz', 'mass',
        'atomic_number', 'formal_charge',
        'name')

class Residue(Handle):
    __slots__ = ()

    def __repr__(self):
        d=self.data()
        return '<Residue %s %d%s>' % (d.name, d.resid, d.insertion)

    def data(self): 
        return self._ptr.residue(self._id)

    def remove(self):
        ''' remove this Residue from the System '''
        self._ptr.delResidue(self._id)

    def addAtom(self):
        ''' append a new Atom to this Residue and return it '''
        return Atom(self._ptr, self._ptr.addAtom(self._id))

    @property
    def atoms(self):
        ''' list of Atoms in this Residue '''
        return [Atom(self._ptr, i) for i in self._ptr.atomsForResidue(self._id)]

    @property
    def natoms(self):
        ''' number of atoms in this residue '''
        return self._ptr.atomCountForResidue(self._id)

    @property
    def chain(self):
        ''' parent chain '''
        return Chain(self._ptr, self.data().chain)
    @chain.setter
    def chain(self, chn): self._ptr.setChain(self._id, chn.id)

    @property
    def resid(self):
        ''' the PDB residue identifier '''
        return self.data().resid
    @resid.setter
    def resid(self, val):
        self.data().resid = val

    @property
    def name(self):
        ''' residue name '''
        return self.data().name
    @name.setter
    def name(self, val):
        self.data().name = val

    @property
    def insertion(self):
        ''' insertion code '''
        return self.data().insertion
    @insertion.setter
    def insertion(self, val):
        self.data().insertion = val

    @property
    def center(self):
        ''' return geometric center of positions of atoms in residue '''
        pos=self._ptr.getPositions(self._ptr.atomsForResidue(self._id))
        return numpy.mean(pos,axis=0)

    def selectAtom(self, name=None):
        ''' Returns a single Atom from this residue with the given name,
        or None if no such atom is present.  If multiple atoms in the
        residue have that name, raise an exception. '''
        atoms = [a for a in self.atoms if (name is None or a.name==name)]
        if not atoms: return None
        if len(atoms) is 1: return atoms[0]
        raise ValueError, "Found %d atoms with given name" % (len(atoms))

class Chain(Handle):
    __slots__ = ()

    def data(self): return self._ptr.chain(self._id)

    def remove(self):
        ''' remove this Chain from the System '''
        self._ptr.delChain(self._id)

    def addResidue(self):
        ''' append a new Residue to this Chain and return it '''
        return Residue(self._ptr, self._ptr.addResidue(self._id))

    @property
    def residues(self):
        ''' list of Residues in this Chain '''
        return [Residue(self._ptr, i) for i in self._ptr.residuesForChain(self._id)]

    @property
    def nresidues(self):
        ''' number of residues in this chain '''
        return self._ptr.residueCountForChain(self._id)

    def selectResidue(self, resid=None, name=None, insertion=None):
        ''' Returns a single Residue with the given resid, name, and/or
        insertion code.  If no such residue is found, returns None.  If
        multiple such residues are found within this chain, raises an
        exception. '''
        residues = [r for r in self.residues if 
                    (resid is None or r.resid==resid) and
                    (name is None or r.name==name) and
                    (insertion is None or r.insertion==insertion)]
        if not residues: return None
        if len(residues) is 1: return residues[0]
        raise ValueError, "Found %d residues with given resid, name or insertion" % (len(residues))

__add_properties(Chain, 'name', 'segid')


class Param(Handle):
    __slots__=()

    def __repr__(self): return '<Param %d>' % self._id

    @property
    def id(self): 
        ''' id in parent table '''
        return self._id
    
    @property
    def table(self): 
        ''' parent ParamTable '''
        return ParamTable(self._ptr)

    def keys(self):
        ''' sorted list of available properties '''
        return sorted(self.table.props)

    def __setitem__(self, prop, val):
        ''' update the value of prop with val '''
        p=self._ptr
        col=p.propIndex(prop)
        if col==_msys.BadId:
            raise KeyError, "No such property '%s'" % prop
        p.setProp(self._id, col, val)

    def __getitem__(self, prop):
        ''' get the value of prop '''
        p=self._ptr
        col=p.propIndex(prop)
        if col==_msys.BadId:
            raise KeyError, "No such property '%s'" % prop
        return p.getProp(self._id, col)

    def duplicate(self):
        ''' create a new entry in the parent parameter table with the
        same values as this one, returning it. '''
        return Param(self._ptr, self._ptr.duplicate(self._id))

class ParamTable(object):
    __slots__=('_ptr',)

    def __init__(self, _ptr):
        self._ptr = _ptr
        ''' Construct from ParamTablePtr.
        Do not invoke directly; use CreateParamTable() instead.  '''

    def __eq__(self, x): 
        return self.__class__==type(x) and self._ptr==x._ptr

    def __ne__(self, x): 
        return self.__class__!=type(x) or  self._ptr!=x._ptr

    def __hash__(self): return self._ptr.__hash__()

    def addParam(self):
        ''' add and return a new Param() '''
        return Param(self._ptr, self._ptr.addParam())

    def addProp(self, name, type):
        ''' add a new property of the given type, which must be int,
        float, or str.  '''
        self._ptr.addProp(name,type)

    @property
    def props(self):
        ''' names of the properties '''
        p=self._ptr
        return [p.propName(i) for i in range(p.propCount())]

    @property
    def nprops(self):
        ''' number of properties '''
        return self._ptr.propCount()

    def propType(self, name):
        ''' type of the property with the given name '''
        p=self._ptr
        return p.propType(p.propIndex(name))

    def delProp(self, name):
        ''' removes the property with the given name. ''' 
        self._ptr.delProp(self._ptr.propIndex(name))

    def param(self, id):
        ''' fetch the Param with the given id '''
        return Param(self._ptr, id)

    @property
    def nparams(self):
        ''' number of Params '''
        return self._ptr.paramCount()

    @property
    def params(self):
        ''' list of all Params in table '''
        return [Param(self._ptr, i) for i in self._ptr.params()]

    def find(self, name, value):
        ''' return the Params with the given value for name '''
        proptype = self.propType(name)
        col = self._ptr.propIndex(name)
        value = proptype(value)
        if proptype is int:
            f=self._ptr.findInt
        elif proptype is float:
            f=self._ptr.findFloat
        else:
            f=self._ptr.findString
        return [self.param(x) for x in f(col, value)]
        
class Term(Handle):
    __slots__=()

    def __repr__(self): return '<Term %d>' % self._id

    @property
    def paramid(self):
        return self._ptr.param(self._id)
    @paramid.setter
    def paramid(self, id):
        self._ptr.setParam(self._id, id)

    def remove(self):
        ''' remove the given Term from its TermTable '''
        self._ptr.delTerm(self._id)

    @property
    def id(self): 
        ''' id of this term in its TermTable '''
        return self._id

    @property
    def param(self): 
        ''' The Param corresponding to this Term's parameters '''
        id=self.paramid
        if _msys.bad(id): return None
        return Param(self._ptr.params(), id)
    @param.setter
    def param(self, val):
        if val is None: 
            id = None
        else: 
            if val.table != self.table.params:
                raise RuntimeError, "param comes from a different ParamTable"
            id = val.id
        self.paramid = id

    @property
    def atoms(self):
        ''' list of Atoms for this Term '''
        return [Atom(self._ptr.system(), i) for i in self._ptr.atoms(self._id)]

    @property
    def table(self): 
        ''' parent TermTable '''
        return TermTable(self._ptr)

    def keys(self):
        ''' union of table.params.props and table.term_props '''
        return sorted(set(self.table.params.props).union(self.table.term_props))

    def __getitem__(self, prop):
        ''' get the value of property prop '''
        ptr = self._ptr
        # first try term properties
        col = ptr.termPropIndex(prop)
        if not _msys.bad(col):
            return ptr.getTermProp(self._id, col)
        # otherwise use param properties
        p = ptr.params()
        col = p.propIndex(prop)
        if col==_msys.BadId:
            raise KeyError, "No such property '%s'" % prop
        id=self.paramid
        if id==_msys.BadId:
            # The user asked for a valid property, but the term doesn't
            # have an assigned param.  
            raise RuntimeError, "No assigned param for '%s'" % repr(self)
        return p.getProp(id, col)

    def __setitem__(self, prop, val):
        ''' set the value of property prop '''
        ptr = self._ptr
        # first try term properties
        col = ptr.termPropIndex(prop)
        if not _msys.bad(col):
            return ptr.setTermProp(self._id, col, val)
        # otherwise use param properties, duplicating if necessary
        p = ptr.params()
        col=p.propIndex(prop)
        if col==_msys.BadId:
            raise KeyError, "No such property '%s'" % prop
        id=self.paramid
        if id==_msys.BadId:
            raise RuntimeError, "No assigned param for '%s'" % repr(self)
        if p.refcount(id) > 1:
            id=p.duplicate(id)
            self.paramid=id
        p.setProp(id, col, val)

class TermTable(object):

    __slots__=('_ptr',)

    def __init__(self, _ptr):
        ''' Construct from TermTablePtr.
        Do not invoke directly; use System.addTable or System.table instead '''
        self._ptr = _ptr

    def __eq__(self, x): return self._ptr==x._ptr
    def __ne__(self, x): return self._ptr!=x._ptr
    def __hash__(self): return self._ptr.__hash__()

    def remove(self):
        ''' Remove this table from its parent system '''
        self._ptr.destroy()

    def coalesce(self):
        ''' Reassign param for each Term in this Table to a member
        of the distinct set of Params used by those Terms.
        '''
        self._ptr.coalesce()

    @property
    def name(self):
        ''' name of this table '''
        return self._ptr.name()
    @name.setter
    def name(self, newname):
        self._ptr.rename(newname)

    def __repr__(self):
        return "<TermTable '%s'>" % self.name

    @property
    def params(self): 
        ''' The ParamTable for terms in this table. '''
        return ParamTable(self._ptr.params())

    @params.setter
    def params(self, p):
        self._ptr.resetParams(p._ptr)

    @property
    def system(self): 
        ''' The System whose atoms are referenced by this table. '''
        sys=self._ptr.system()
        if sys is None: return None
        return System(sys)

    @property
    def term_props(self): 
        ''' names of the custom properties '''
        p=self._ptr
        return [p.termPropName(i) for i in range(p.termPropCount())]

    @property
    def natoms(self): 
        ''' number of atoms in each term '''
        return self._ptr.atomCount()

    @property
    def category(self): 
        ''' A string describing what kind of TermTable this is.  
        Possibilities are: *bond*, *constraint*, *virtual*, *polar*, 
        *nonbonded*, and *exclusion*.'''
        return _msys.print_category(self._ptr.category)
    @category.setter
    def category(self, val): 
        self._ptr.category=_msys.parse_category(val)

    @property
    def nterms(self): 
        ''' number of terms '''
        return self._ptr.termCount()

    def delTermsWithAtom(self, atom):
        ''' remove all terms whose atoms list contains the given Atom '''
        assert atom.system == self.system, "atom is from different System"
        self._ptr.delTermsWithAtom(atom.id)

    def findWithAll(self, atoms):
        ''' return the terms that contain all the given atoms in any order '''
        if not atoms: return []
        ptr, ids = _find_ids(atoms)
        if ptr!=self.system._ptr:
            raise ValueError, "atoms are from a different System"
        return [self.term(x) for x in self._ptr.findWithAll(ids)]

    def findWithAny(self, atoms):
        ''' return the terms that contain at least one of the given atoms '''
        if not atoms: return []
        ptr, ids = _find_ids(atoms)
        if ptr!=self.system._ptr:
            raise ValueError, "atoms are from a different System"
        return [self.term(x) for x in self._ptr.findWithAny(ids)]

    def findWithOnly(self, atoms):
        ''' return the terms that contain only the given atoms '''
        if not atoms: return []
        ptr, ids = _find_ids(atoms)
        if ptr!=self.system._ptr:
            raise ValueError, "atoms are from a different System"
        return [self.term(x) for x in self._ptr.findWithOnly(ids)]

    def findExact(self, atoms):
        ''' return the terms that contain precisely the given atoms in the 
        given order. '''
        if not atoms: return []
        ptr, ids = _find_ids(atoms)
        if ptr!=self.system._ptr:
            raise ValueError, "atoms are from a different System"
        return [self.term(x) for x in self._ptr.findExact(ids)]

    @property
    def terms(self):
        ''' returns a list of all the Terms in the table '''
        return [Term(self._ptr, i) for i in self._ptr.terms()]

    def term(self, id):
        ''' returns the Term in the table with the given id '''
        return Term(self._ptr, id)

    def hasTerm(self, id):
        ''' Does a Term with the given id exist in the table? '''
        return self._ptr.hasTerm(id)

    def addTermProp(self, name, type):
        ''' add a custom Term property of the given type '''
        self._ptr.addTermProp(name, type)

    def delTermProp(self, name):
        ''' remove the custom Term property '''
        self._ptr.delTermProp(self._ptr.termPropIndex(name))

    def termPropType(self, name):
        ''' type of the given Term property '''
        p=self._ptr
        return p.termPropType(p.termPropIndex(name))

    def addTerm(self, atoms, param = None):
        ''' Add a Term to the table, with given initial param.  The atoms
        list must have natoms elements, and each Atom must belong to the 
        same System as the TermTable.  If param is not None, it must
        belong to the ParamTable held by the TermTable.
        '''
        if param is not None:
            if param.table != self.params:
                raise RuntimeError, "param comes from a different ParamTable"
            param = param.id
        n=len(atoms)
        ids=[None]*n
        for i in range(n):
            a=atoms[i]
            ids[i]=a.id
            if a.system != self.system:
                raise RuntimeError, "Cannot add atoms from different system"
        return Term(self._ptr, self._ptr.addTerm(ids, param))

    @property
    def override_params(self):
        ''' parameter table containing override values '''
        return ParamTable(self._ptr.overrides().params())

    @property
    def noverrides(self):
        ''' number of parameter overrides '''
        return self._ptr.overrides().count()

    def overrides(self):
        ''' return a mapping from pairs of Params in self.params to a Param
        in self.override_params. '''
        d=dict()
        p = self.params
        o = self._ptr.overrides()
        op = self.override_params
        for i,j in o.list():
            pi = p.param(i)
            pj = p.param(j)
            d[(pi,pj)] = op.param(o.get(i,j))
        return d

    def setOverride(self, pi, pj, op):
        ''' override the interaction between params pi and pj with the
        interaction given by op.  pi and pj must be Params from self.params;
        op must be a param from self.override_params, or None to remove
        the override.
        '''
        if pi.__class__ is not Param: raise TypeError, "pi must be Param"
        if pj.__class__ is not Param: raise TypeError, "pj must be Param"
        if pi.table != self.params: 
            raise ValueError, "pi must be from self.params"
        if pj.table != self.params: 
            raise ValueError, "pj must be from self.params"
        if op is None:
            return self._ptr.overrides().del_(pi.id, pj.id)
        if op.__class__ is not Param: 
            raise TypeError, "op must be Param or None"
        if op.table != self.override_params: 
            raise ValueError, "op must be from self.override_params"
        self._ptr.overrides().set(pi.id, pj.id, op.id)

    def getOverride(self, pi, pj):
        ''' get override for given pair of params, or None if not present. '''
        if pi.__class__ is not Param: raise TypeError, "pi must be Param"
        if pj.__class__ is not Param: raise TypeError, "pj must be Param"
        if pi.table != self.params: 
            raise ValueError, "pi must be from self.params"
        if pj.table != self.params: 
            raise ValueError, "pj must be from self.params"
        id = self._ptr.overrides().get(pi.id, pj.id)
        if _msys.bad(id): return None
        return self.override_params.param(id)

class System(object):

    __slots__ = ('_ptr', '_atoms')

    def __init__(self, _ptr):
        ''' Construct from SystemPtr.
        Do not invoke directly; use CreateSystem() instead.
        '''
        self._ptr = _ptr
        self._atoms = []

    def __eq__(self, x): 
        return self.__class__==type(x) and self._ptr==x._ptr

    def __ne__(self, x): 
        return self.__class__!=type(x) or  self._ptr!=x._ptr

    def __hash__(self): return self._ptr.__hash__()

    def __repr__(self): return "<System '%s'>" % self.name

    @property
    def name(self): 
        ''' The name of the System, taken from the input file name '''
        return self._ptr.name
    @name.setter
    def name(self, s): self._ptr.name=s

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
        ''' remove all Atoms in elems from the System '''
        self._ptr.delAtoms(e.id for e in elems)

    def delBonds(self, elems):
        ''' remove all Bonds in elems from the System '''
        self._ptr.delBonds(e.id for e in elems)

    def delResidues(self, elems):
        ''' remove all Residues in elems from the System '''
        self._ptr.delResidues(e.id for e in elems)

    def delChains(self, elems):
        ''' remove all Chains in elems from the System '''
        self._ptr.delChains(e.id for e in elems)

    def atom(self, id):
        ''' return the atom with the specified id '''
        return Atom(self._ptr, id)

    def bond(self, id):
        ''' return the bond with the specified id '''
        return Bond(self._ptr, id)

    def residue(self, id):
        ''' return the residue with the specified id '''
        return Residue(self._ptr, id)

    def chain(self, id):
        ''' return the chain with the specified id '''
        return Chain(self._ptr, id)

    def findBond(self, a1, a2):
        ''' return the bond between the specified atoms, or None if not found
        '''
        assert a1.system == self
        assert a2.system == self
        id=self._ptr.findBond(a1.id, a2.id)
        if _msys.bad(id): return None
        return Bond(self._ptr, id)

    @property
    def cell(self):
        ''' The GlobalCell for this System '''
        return self._ptr.global_cell

    @property
    def nonbonded_info(self):
        ''' NonbondedInfo for this System '''
        return self._ptr.nonbonded_info

    @property
    def natoms(self):
        ''' number of atoms '''
        return self._ptr.atomCount()

    @property
    def nbonds(self):
        ''' number of bonds '''
        return self._ptr.bondCount()

    @property
    def nresidues(self):
        ''' number of residues '''
        return self._ptr.residueCount()

    @property
    def nchains(self):
        ''' number of chains '''
        return self._ptr.chainCount()

    @property
    def atoms(self):
        ''' return list of all atoms in the system '''
        atms=self._update_atoms()
        return [atms[i] for i in self._ptr.atomsAsList()]

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
        self._ptr.addAtomProp(name, type)

    def delAtomProp(self, name):
        ''' remove the given custom atom property '''
        self._ptr.delAtomProp(self._ptr.atomPropIndex(name))


    @property
    def atom_props(self):
        ''' return the list of custom atom properties. '''
        p=self._ptr
        return [p.atomPropName(i) for i in range(p.atomPropCount())]

    def atomPropType(self, name):
        ''' type of the given atom property '''
        return self._ptr.atomPropType(self._ptr.atomPropIndex(name))

    def atomsGroupedBy(self, prop):
        ''' Return dictionary mapping representative values of the given
        atom property to lists of atoms having that property.  If the
        property does not exist in this system, returns an empty dictionary.
        '''
        d=dict()
        if hasattr(Atom, prop):
            getter=lambda x: getattr(x,prop)
        elif prop in self.atom_props:
            getter=lambda x: x[prop]
        else:
            return d
        atms=self._update_atoms()
        for a in atms:
            key = getter(a)
            d.setdefault(key,[]).append(a)
        return d

    @property
    def positions(self):
        ''' Nx3 list of lists of positions of all atoms '''
        return self._ptr.getPositions()
    @positions.setter
    def positions(self, pos):
        self._ptr.setPositions(pos)

    def getPositions(self):
        ''' get copy of positions as Nx3 array '''
        return self._ptr.getPositions()

    def setPositions(self, pos):
        ''' set positions from Nx3 array '''
        self._ptr.setPositions(pos)

    def getVelocities(self):
        ''' get copy of velocities as N3x array '''
        return self._ptr.getVelocities()

    def setVelocities(self, vel):
        ''' set velocities from Nx3 array '''
        self._ptr.setVelocities(vel)

    def setCell(self, cell):
        ''' set unit cell from from 3x3 array '''
        for i in range(3):
            self.cell[i][:] = map(float, cell[i])

    def getCell(self):
        ''' return unit cell as 3x3 array '''
        return [self.cell[i][:] for i in range(3)]

    @property
    def center(self):
        ''' return geometric center of positions of all atoms '''
        return numpy.mean(self.getPositions(), axis=0)

    def translate(self, xyz):
        ''' shift coordinates by given amount '''
        x,y,z = xyz
        for id in self._ptr.atoms():
            atm=self._ptr.atom(id)
            atm.x += x
            atm.y += y
            atm.z += z

    ###
    ### bond properties
    ###

    def addBondProp(self, name, type):
        ''' add a custom bond property with the given name and type.
        type should be int, float, or str.
        '''
        self._ptr.addBondProp(name, type)

    def delBondProp(self, name):
        ''' remove the given custom bond property '''
        self._ptr.delBondProp(self._ptr.bondPropIndex(name))


    @property
    def bond_props(self):
        ''' return the list of custom bond properties. '''
        p=self._ptr
        return [p.bondPropName(i) for i in range(p.bondPropCount())]

    def bondPropType(self, name):
        ''' type of the given bond property '''
        return self._ptr.bondPropType(self._ptr.bondPropIndex(name))

    
    @property
    def topology(self, with_glue=True):
        ''' list of bonded atoms for each atom in the System, augmented
        by whatever glue may be present.
        '''
        top = [[b.id for b in a.bonded_atoms] for a in self.atoms]
        if with_glue and 'glue' in self.auxtable_names:
            gtable = self.auxtable('glue')
            for t in gtable.params:
                p0 = t['p0']
                p1 = t['p1']
                if p1 not in top[p0]:
                    top[p0].append(p1)
                if p0 not in top[p1]:
                    top[p1].append(p0)
        return top

    ###
    ### operations on term tables
    ###

    @property
    def table_names(self):
        ''' names of the tables in the System '''
        return [x for x in self._ptr.tableNames()]

    @property
    def tables(self):
        ''' all the tables in the System '''
        return [TermTable(self._ptr.table(x)) for x in self._ptr.tableNames()]

    def table(self, name):
        ''' Get the TermTable with the given name '''
        ptr=self._ptr.table(name)
        if ptr is None:
            raise ValueError, "No such table '%s'" % name
        return TermTable(ptr)

    def addTable(self, name, natoms, params = None):
        ''' add a table with the given name and number of atoms.
        If a table with the same name already exists, it is returned,
        otherwise the newly created table is returned.  If no ParamTable
        params is supplied, a new one is created.  '''
        if params is None: params = CreateParamTable()
        name = str(name)
        natoms = int(natoms)
        return TermTable(self._ptr.addTable(name, natoms, params._ptr))

    def addTableFromSchema(self, type, name = None):
        ''' Add a table to the system if it not already present, 
        returning it.  If optional name field is provided, the table
        will be added with the given name; otherwise the name is taken
        from the table schema. '''
        if name is None: name=type
        return TermTable(self._ptr.addTableFromSchema(type,name))

    def coalesceTables(self):
        ''' Invoke TermTable.coalesce on each table '''
        self._ptr.coalesceTables()


    ###
    ### auxiliary tables
    ###

    @property
    def auxtable_names(self):
        ''' names of the auxiliary tables '''
        return [x for x in self._ptr.auxTableNames()]

    @property
    def auxtables(self):
        ''' all the auxiliary tables '''
        return [ParamTable(self._ptr.auxTable(x)) for x in self._ptr.auxTableNames()]

    def auxtable(self, name):
        ''' auxiliary table with the given name '''
        ptr=self._ptr.auxTable(name)
        if ptr is None:
            raise ValueError, "No such extra table '%s'" % name
        return ParamTable(ptr)

    def addAuxTable(self, name, table):
        ''' add or replace extra table with the given name. '''
        self._ptr.addAuxTable(name, table._ptr)

    def delAuxTable(self, name):
        ''' remove auxiliary table with the given name. '''
        self._ptr.delAuxTable(name)

    def addNonbondedFromSchema(self, funct, rule=""):
        '''
        Add a nonbonded table to the system, and configure the nonbonded
        info according to funct and rule.  funct must be the name of recognized
        nonbonded type.  rule is not checked; at some point in the future we
        might start requiring that it be one of the valid combining rules for
        the specified funct.  If nonbonded_info's vdw_funct and vdw_rule
        are empty, they are overridden by the provided values; otherwise, the
        corresponding values must agree if funct and rule are not empty.
        A nonbonded table is returned.
        '''
        return TermTable(self._ptr.addNonbondedFromSchema(funct,rule))

    def _update_atoms(self):
        p=self._ptr
        atms=self._atoms
        n=len(atms)
        A=Atom
        for i in xrange(n,p.maxAtomId()):
            atms.append(A(p,i))
        return atms

    def atomsel(self, seltext):
        ''' Create and return an atom selection object (Atomsel). '''
        return Atomsel(self._ptr, self._ptr.select(seltext))

    def select(self, seltext):
        ''' return a list of Atoms satisfying the given VMD atom selection. '''
        ids=self._ptr.selectAsList(seltext)
        atms=self._update_atoms()
        return [atms[i] for i in ids]

    def selectIds(self, seltext):
        ''' Return the ids of the Atoms satisfying the given VMD atom
        selection.  This can be considerably faster than calling select().
        '''
        return self._ptr.selectAsList(seltext)

    def selectChain(self, name=None, segid=None):
        ''' Returns a single Chain with the matching name and/or segid,
        or raises an exception if no single such chain is present.
        '''
        chains = [c for c in self.chains if 
                    (name is None or c.name==name) and
                    (segid is None or c.segid==segid)]
        if not chains: return None
        if len(chains) is 1: return chains[0]
        raise ValueError, "Found %d chains with given name and segid" % (len(chains))

    def append(self, system):
        ''' Appends atoms and forcefield from system to self.  Returns
        a list of of the new created atoms in self.  Systems must have
        identical nonbonded_info.vdw_funct. '''
        p=self._ptr
        ids=p.append(system._ptr)
        return [Atom(p,i) for i in ids]

    def clone(self, seltext=None):
        ''' Clone the System, returning a new System.  If seltext is provided,
        it should be a valid VMD atom selection, and only the selected atoms 
        will be cloned.

        This operation is equivalent to "msys.CloneSystem(self.select(seltext))"
        where seltext defaults to 'all'.
        '''
        ptr = self._ptr
        if seltext is None:
            ids = ptr.atoms()
        else:
            ids = ptr.select(seltext)
        return System( _msys.Clone(ptr, ids))

    def sorted(self):
        ''' Return a clone of the system with atoms reordered based on their 
        order of appearance in a depth-first traversal of the structure 
        hierarchy.
        '''
        ptr=self._ptr
        return System(_msys.Clone(ptr, ptr.orderedIds()))

    def permuted(self, perm):
        ''' Return a permutation of the atoms in the system.  perm should
        be either a permutation of range(natoms), or a permutation of the
        atoms in the system.
        '''
        if len(perm)!=self.natoms:
            raise ValueError, "perm has %d elements, expected natoms=%d" % (
                    len(perm), self.natoms)
        if isinstance(perm[0], Atom):
            ptr, ids = _find_ids(perm)
            if ptr != self._ptr:
                raise ValueError, "Atoms in perm are not from this System"
        else:
            ptr = self._ptr
            ids = _msys.IdList()
            for p in perm: ids.append(perm.id)
        if not ptr.validPermutation(ids):
            raise ValueError, "Not a valid permutation"
        return System(_msys.Clone(self._ptr, ids))


    def guessBonds(self, replace=True, reanalyze=True):
        ''' Guess bond connectivity based on an atomic-number based
        atom radius.  

        Replaces any existing bonds, unless replace=False is specified.

        Reanalyzes fragids and atom types unless reanalyze=False is specified.
        In that case, you MUST call updateFragids() manually before making
        any use of the fragment assignment (fragids will be out of date).
        '''
        ptr=self._ptr
        if replace:
            ptr.delBonds(ptr.bonds())
        _msys.GuessBondConnectivity(ptr)
        if reanalyze:
            ptr.analyze()

    def analyze(self):
        ''' Assign atom and residue types.  This needs to be called
        manually only if you create a system from scratch, using 
        msys.CreateSystem(); in that case, analyze() should be called
        before performing any atom selections.
        '''
        self._ptr.analyze()

    def updateFragids(self):
        ''' Find connected sets of atoms, and assign each a 0-based id,
        stored in the fragment property of the atom.  Return a list of
        fragments as a list of lists of atoms. '''
        p = self._ptr
        frags = p.updateFragids()
        result=[]
        for f in frags:
            result.append( [Atom(p,i) for i in f] )
        return result

    @property
    def provenance(self):
        ''' return a list of Provenance entries for this system '''
        return self._ptr.provenance()

    def addSelectionMacro(self, macro, definition):
        ''' Add and/or replace the given atom selection macro with the given
        definition, which must be parseable.
        '''
        self._ptr.addSelectionMacro(str(macro).strip(), str(definition))

    def delSelectionMacro(self, macro):
        ''' Remove the given macro from the dictionary for this system '''
        self._ptr.delSelectionMacro(str(macro))

    @property
    def selection_macros(self):
        ''' Return a list of the selection macros defined for this system.
        '''
        return [x for x in self._ptr.selectionMacros()]

    def selectionMacroDefinition(self, macro):
        ''' Return the definition for the given macro, or None if not defined
        '''
        x = self._ptr.selectionMacroDefinition(str(macro).strip())
        return x if x else None

    def gluePairs(self):
        ''' return a list of the glue pairs as tuples '''
        return self._ptr.gluePairs()

    def addGluePair(self, p0, p1):
        ''' Add a pair of atoms to the glue.  '''
        return self._ptr.addGluePair(p0,p1)

    def hasGluePair(self, p0, p1):
        ''' Does the pair of atoms exist in the glue? '''
        return self._ptr.hasGluePair(p0,p1)

    def delGluePair(self, p0, p1):
        ''' Remove a pair of atoms from the glue.  '''
        return self._ptr.delGluePair(p0,p1)

class SmartsPattern(object):
    ''' A class representing a compiled SMARTS pattern '''
    def __init__(self, pattern):
        ''' Initialize with SMARTS pattern '''
        self._pat = _msys.SmartsPattern(str(pattern))

    @property
    def natoms(self):
        ''' Number of atoms in the compiled smarts pattern '''
        return self._pat.atomCount()

    @property
    def pattern(self):
        ''' The pattern used to initialize the object '''
        return self._pat.pattern()

    @property
    def warnings(self):
        ''' Warnings, if any, emitted during compilation '''
        return self._pat.warnings()

    def __repr__(self):
        return "<SmartsPattern '%s'>" % self.pattern

    def findMatches(self, system_or_atoms):
        ''' return list of lists representing ids of matches of this
        pattern that start with the given set of atoms.  Formal charges
        and bond orders must already be present in the system; use
        AssignBondOrderAndFormalCharge to let msys do the assignment.  '''
        if isinstance(system_or_atoms, System):
            ptr = system_or_atoms._ptr
            ids = ptr.atoms()
        else:
            ptr, ids = _find_ids(system_or_atoms)
        return self._pat.findMatches(ptr, ids)

def CreateSystem():
    ''' Create a new, empty System '''
    return System(_msys.SystemPtr.create())

def _find_ids(atoms):
    ''' return the SystemPtr and IdList for the given set of atoms.
    Raise ValueError if atoms come from multiple systems.  Return
    None, [] if the input list is empty.
    '''
    if not atoms:
        return None, []
    ids = _msys.IdList()
    sys = set()
    for a in atoms:
        ids.append(a.id)
        sys.add(a.system)
    if len(sys) > 1:
        raise ValueError, "Input atoms come from multiple systems"
    ptr = sys.pop()._ptr
    return ptr, ids

def CloneSystem(atoms):
    ''' create a new system from the given list of atoms, which must all
    be from the same System and have different ids.  The order of the atoms
    in the new system will be that of the input atoms.  '''
    ptr, ids = _find_ids(atoms)
    if ptr is None:
        raise ValueError, "Cannot clone an empty list of atoms"
    return System( _msys.Clone(ptr, ids))

def CreateParamTable():
    ''' Create a new, empty ParamTable '''
    return ParamTable(_msys.ParamTablePtr.create())

def LoadDMS(path=None, structure_only=False, buffer=None ):
    ''' Load the DMS file at the given path and return a System containing it.
    If structure_only is True, only Atoms, Bonds, Residues and Chains will
    be loaded, along with the GlobalCell, and no pseudos (atoms with atomic 
    number less than one) will be loaded.

    If the buffer argument is provided, it is expected to hold the contents
    of a DMS file, and the path argument will be ignored.
    '''
    if buffer is None and path is None:
        raise ValueError, "Must provide either path or buffer"
    if buffer is not None and path is not None:
        raise ValueError, "Must provide either path or buffer"

    if path is not None:
       ptr = _msys.ImportDMS(path, structure_only )
    else:
       ptr = _msys.ImportDMSFromBuffer(buffer, structure_only)
    return System(ptr)


def LoadMAE(path=None, ignore_unrecognized = False, buffer=None,
                       structure_only = False):
    ''' load the MAE file at the given path and return a System containing it.
    Forcefield tables will be created that attempt to match as closely as
    possible the force terms in the MAE file; numerical differences are bound
    to exist.  If ignore_unrecognized is True, ignore unrecognized force
    tables.

    If the buffer argument is provided, it is expected to hold the contents
    of a DMS file, and the path argument will be ignored.

    If the contents of the file specified by path, or the contents of buffer,
    are recognized as being gzip-compressed, they will be decompressed on
    the fly.
    
    If structure_only is True, no forcefield components will be loaded.  '''

    if buffer is None and path is None:
        raise ValueError, "Must provide either path or buffer"
    if buffer is not None and path is not None:
        raise ValueError, "Must provide either path or buffer"
    ignore_unrecognized = bool(ignore_unrecognized)
    structure_only = bool(structure_only)

    if path is not None:
        ptr = _msys.ImportMAE(path, ignore_unrecognized, structure_only )
    else:
        ptr = _msys.ImportMAEFromBuffer( buffer, ignore_unrecognized, 
                                                 structure_only )
    return System(ptr)

def LoadPDB(path):
    ''' Load a PDB file at the given path and return a System.
    No bonds will be created, even if CONECT records are parent.
    '''
    path=str(path)
    ptr = _msys.ImportPDB(path)
    return System(ptr)

def LoadPrmTop(path, structure_only = False):
    ''' Load an Amber7 prmtop file at the given path and return a System.
    Coordinates and global cell information are not present in the file.
    '''
    path=str(path)
    ptr = _msys.ImportPrmTop(path, bool(structure_only))
    return System(ptr)

def LoadMol2(path, multiple=False):
    ''' Load a mol2 file at the given path. If multiple is True, return a list
    of Systems, one for each MOLECULE record.  Otherwise, return just one
    System corresponding to the first MOLECULE record.
    '''
    if multiple:
        return [System(p) for p in _msys.ImportMOL2Many(path)]
    return System(_msys.ImportMOL2(path))

def LoadXYZ(path):
    ''' Load an xyz file at the given path.  Guesses bonds based on 
    guessed atomic numbers based on atom name.
    '''
    return System(_msys.ImportXYZ(path))


def Load(path, structure_only = False):
    ''' Infer the file type of path and load the file.
    Returns a new System.
    '''
    ptr = _msys.Load(path, bool(structure_only))
    if not ptr:
        raise ValueError, "Could not guess file type of '%s'" % path
    return System(ptr)

def LoadMany(path, structure_only=False):
    ''' Iterate over structures in a file, if the file type supports
    iteration.  

    for mol in LoadMany('input.mol2'): ...
    '''
    it = _msys.LoadIterator.create(str(path), bool(structure_only))
    while True:
        mol = it.next()
        if mol is None:
            break
        yield System(mol)

def ReadCrdCoordinates(mol, path):
    ''' Read coordinates from the given Amber crd file into the given 
    System. 
    '''
    path=str(path)
    if not isinstance(mol, System): raise TypeError, "mol must be a System"
    _msys.ImportCrdCoordinates(mol._ptr, path)

def SaveDMS(system, path):
    ''' Export the System to a DMS file at the given path. '''
    _msys.ExportDMS(system._ptr, path, _msys.Provenance.fromArgs(sys.argv))


def SaveMAE(system, path, with_forcefield = True, with_compression=False ):
    ''' Export the System to an MAE file at the given path. 
    If system is a list, then each element should be a system, and a multi-
    ct MAE file will be written.  '''
    path = str(path)
    prov = _msys.Provenance.fromArgs(sys.argv)
    ff = bool(with_forcefield)
    compress = bool(with_compression)
    try:
        cts = list(system)
    except TypeError:
        _msys.ExportMAE(system._ptr, path, prov, ff, compress) 
    else:
        _msys.ExportMAEMany([x._ptr for x in cts], path, prov, ff, compress)

def SavePDB(system, path):
    ''' Export the System to a PDB file at the given path. '''
    _msys.ExportPDB(system._ptr, str(path))

def SaveMol2(system, path):
    ''' Export the System to a mol2 file at the given path.  You should
    probably call AssignBondOrderAndFormalCharge() before exporting
    the system.  '''
    ptr=system._ptr
    path=str(path)
    prov=_msys.Provenance.fromArgs(sys.argv)
    _msys.ExportMOL2(ptr,path,prov)

def SaveSDF(system, path_or_file):
    ''' Export the system to the given path or file.  If path_or_file is
    a string, it is assumed to be a path; if None, the bytes are returned;
    otherwise, it is assumed to be a file object implementing the write()
    method.
    '''
    s = _msys.ExportSDFBytes(system._ptr)
    if path_or_file is None:
        return s
    if isinstance(path_or_file, str):
        f = file(path_or_file, 'w')
    else:
        f = path_or_file
    f.write(s)
    f.flush()

def TableSchemas():
    ''' available schemas for System.addTableFromSchema '''
    return [s for s in _msys.TableSchemas()]

def NonbondedSchemas():
    ''' available nonbonded schemas for System.addNonbondedFromSchema '''
    return [s for s in _msys.NonbondedSchemas()]

def GetSSSR(atoms, all_relevant=False):
    """Get smallest set of smallest rings (SSSR) for a system fragment.

    The SSSR is in general not unique; the SSSR of a tetrahedron is any
    three of its four triangular faces. The set of rings that is the
    union of all SSSR's (all relevant rings) may be obtained by setting
    all_relevant to True.

    Arguments:
    atoms -- [msys.Atom, ..., msys.Atom] from a single system
    all_relevant -- bool
    Returns: [[msys.Atom, ..., msys.Atom], ..., [msys.Atom, ..., msys.Atom]]
    """
    ptr, ids = _find_ids(atoms)
    rings = _msys.GetSSSR(ptr, ids, all_relevant)
    return [[Atom(ptr, id) for id in ring] for ring in rings]

def AssignSybylTypes(system):
    ''' Assign Sybyl atom and bond types to the given system.  
    Types will be stored in the "sybyl_type" property of each atom and bond.
    Invoke AssignBondOrderAndFormalCharge first if bond orders need to
    be recalculated. '''
    _msys.AssignSybylTypes(system._ptr)

def AssignBondOrderAndFormalCharge(system_or_atoms, total_charge = None):
    """Assign bond orders and formal charges to a molecular system.

    Determines bond orders and formal charges by preferring neutral
    charges and placing negative charges with more electronegative
    atoms, under octet constraints and the total system charge
    constraint. Assigns the bond orders and formal charges to the
    system. Also determines resonance forms and creates and populates a
    'resonant_order' bond property and 'resonant_charge' atom property
    in the system, where resonant order/charge is the average of the
    bond order/formal charge over all resonant forms. Can assign to a
    subset of atoms of the system, provided these atoms form complete
    connected fragments.

    WARNING: calling this function on a chemically incomplete system,
    i.e. just protein backbone, may cause msys to hang indefinitely.

    Arguments:
    system_or_atoms: either a System or a list of Atoms
    total_charge: if not None, integral total charge
    """
    if isinstance(system_or_atoms, System):
        ptr = system_or_atoms._ptr
        ids = ptr.atoms()
    else:
        ptr, ids = _find_ids(system_or_atoms)

    if total_charge is None:
        _msys.AssignBondOrderAndFormalCharge(ptr, ids)
    else:
        _msys.AssignBondOrderAndFormalCharge(ptr, ids, int(total_charge))

class Graph(object):
    @classmethod
    def _from_boost(cls, _ptr, _sys):
        graph = cls(CreateSystem())
        graph._ptr = _ptr
        graph._sys = _sys
        return graph

    def __init__(self, system_or_atoms):
        if isinstance(system_or_atoms, System):
            ptr = system_or_atoms._ptr
            ids = ptr.atoms()
        else:
            ptr, ids = _find_ids(system_or_atoms)
        self._ptr = _msys.GraphPtr.create(ptr,ids)
        self._sys = ptr

    def size(self):
        ''' number of atoms in graph '''
        return self._ptr.size()

    def atoms(self):
        ''' ordered atoms in graph '''
        return [Atom(self._sys, i) for i in self._ptr.atoms()]

    def hash(self):
        ''' string hash of atoms and bonds in graph '''
        return self._ptr.hash()

    def match(self, graph):
        ''' Find a graph isomorphism between self and the given Graph.
        If no isomorphism could be found, return None; otherwise return
        mapping from atoms in this graph to atoms in that graph.
        '''
        if not isinstance(graph, Graph):
            raise TypeError, "graph argument must be an instance of msys.Graph"
        t = self._ptr.match(graph._ptr)
        if t is not None:
            t=dict((Atom(self._sys, i), Atom(graph._sys,j)) for i,j in t)
        return t

    def matchAll(self, graph, substructure=False):
        ''' Find all graph isomorphisms between self and the given Graph.
        If no isomorphism could be found, return empty list; otherwise return
        list of dicts mapping atoms in this graph to atoms in that graph. If
        substructure is True, return isomorphisms between self and any subgraph
        of the given Graph.
        '''
        if not isinstance(graph, Graph):
            raise TypeError, "graph argument must be an instance of msys.Graph"
        t = self._ptr.matchAll(graph._ptr, substructure)
        return [dict((Atom(self._sys, i), Atom(graph._sys, j)) for i,j in item)
                for item in t]

def FindDistinctFragments(system):
    ''' Return fragids of representative fragments.  '''
    return _msys.FindDistinctFragments(system._ptr)

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
Vec3.__repr__=lambda self: str((self.x, self.y, self.z))

''' customize GlobalCell '''
def __globalcell_str(self):
    A=[x for x in self.A]
    B=[x for x in self.B]
    C=[x for x in self.C]
    return str([A,B,C])
GlobalCell.__str__ = __globalcell_str




