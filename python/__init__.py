
'''
This is the high-level Python interface for msys, intended for use
by chemists.
'''

import _msys
import numpy
import sys
import tempfile

from _msys import NonbondedInfo, version, hexversion
from _msys import RadiusForElement, MassForElement, ElementForAbbreviation
from _msys import GuessAtomicNumber, AbbreviationForElement
from _msys import ElectronegativityForElement
from _msys import PeriodForElement, GroupForElement
from _msys import HydrogenBond, FetchPDB
from _msys import BadId
from atomsel import Atomsel
import molfile

class Handle(object):
    __slots__ = ('_ptr', '_id')

    def __init__(self, _ptr, _id):
        self._ptr = _ptr
        self._id = int(_id)

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
    ''' Represents a bond in a System '''

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

    def other(self, atom):
        ''' atom in bond not the same as given atom '''
        return Atom(self._ptr, self.data().other(atom.id))

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

    @property
    def order(self):
        ''' bond order (int) '''
        return self.data().order
    @order.setter
    def order(self, val):
        self.data().order = val

class Atom(Handle):
    ''' Represents an atom (or pseudoparticle) in a chemical system '''
    __slots__ = ()

    def __repr__(self): return '<Atom %d>' % self._id

    @property
    def fullname(self):
        res = self.residue
        chn = res.chain
        return '%s:%s%d:%s' % (chn.name, res.name, res.resid, self.name)

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

    @property
    def nhydrogens(self):
        ''' number of bonded hydrogens '''
        return len([a for a in self.bonded_atoms if a.atomic_number==1])

    @property
    def valence(self):
        ''' sum of bond orders '''
        return sum([b.order for b in self.bonds if b.other(self).atomic_number > 0])

__add_properties(Atom, 
        'fragid', 'x', 'y', 'z', 'charge',
        'vx', 'vy', 'vz', 'mass',
        'atomic_number', 'formal_charge',
        'name', 'aromatic')

class Residue(Handle):
    ''' Represents a residue (group of Atoms) in a System '''
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
    ''' Represents a chain (of Residues) in a System '''
    __slots__ = ()

    def __repr__(self): return '<Chain %s>' % self.name

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

    @property
    def ct(self):
        ''' Return the Ct for this chain '''
        return Ct(self._ptr, self.data().ct)
    @ct.setter
    def ct(self, ct): self._ptr.setCt(self._id, ct.id)

__add_properties(Chain, 'name', 'segid')

class Ct(Handle):
    ''' Represents a list of Chains in a System

    The Ct class exists mainly to provide a separate namespace for chains.
    If you merge two systems each of which has a chain A, you probably
    want the chains to remain separate.  Cts accomplish this.

    The Ct class also provides a key-value namespace for assigning 
    arbitrary properties to Systems.
    '''

    __slots__ = ()
    def data(self): return self._ptr.ct(self._id)

    def remove(self):
        ''' remove this Ct from the System '''
        self._ptr.delCt(self._id)

    def addChain(self):
        ''' append a new Chain to this Ct and return it '''
        return Chain(self._ptr, self._ptr.addChain(self._id))

    @property
    def chains(self):
        ''' list of Chains in this Ct '''
        return [Chain(self._ptr, i) for i in self._ptr.chainsForCt(self._id)]

    @property
    def nchains(self):
        ''' number of Chains in this Ct '''
        return self._ptr.chainCountForCt(self._id)

    @property
    def natoms(self):
        ''' number of Atoms in the Ct '''
        return self._ptr.atomCountForCt(self._id)

    @property
    def atoms(self):
        ''' list of Atoms in this Ct '''
        ptr = self._ptr
        return [Atom(ptr, i) for i in ptr.atomsForCt(self._id)]

    @property
    def bonds(self):
        ''' list of Bonds in this Ct '''
        ptr = self._ptr
        return [Bond(ptr, i) for i in ptr.bondsForCt(self._id)]

    def append(self, system):
        ''' Appends atoms and forcefield from system to self.  Returns
        a list of of the new created atoms in self.  Systems must have
        identical nonbonded_info.vdw_funct.  Does not overwrite the
        global cell information in self.  '''
        p=self._ptr
        ids=p.append(system._ptr, self.id)
        return [Atom(p, i) for i in ids]

    ### TODO: selectChain()

    @property
    def name(self):
        ''' Name of Ct '''
        return self.data().name
    @name.setter
    def name(self, s):
        self.data().name = s

    def keys(self):
        ''' available Ct properties '''
        return self.data().keys()

    def __setitem__(self, key, val):
        ''' set ct property key to val '''
        t = self.data().type(key)
        if t is None:
            t = type(val)
            self.data().add(key, t)
        self.data().set(str(key), t(val))

    def __getitem__(self, key):
        ''' get ct property key '''
        return self.data().get(str(key))

    def __delitem__(self, key):
        ''' remove property key '''
        self.data().remove(str(key))

    def get(self, key, d=None):
        ''' get ct property key, else d, which defaults to None '''
        try:
            return self.data().get(str(key))
        except KeyError:
            return d

class PropertyMap(object):
    def __init__(self, ptr):
        self._ptr = ptr

    def keys(self):
        return self._ptr.keys()

    def __getitem__(self, key):
        return self._ptr.get(str(key))

    def __setitem__(self, key, val):
        self._ptr.set(str(key), type(val), val)

    def __delitem__(self, key):
        self._ptr._del(str(key))


class Param(Handle):
    '''
    A `Param` instance is a reference to a row in a `ParamTable`.  Use the
    ``dict``-style interface to get and set values in the row.  Msys will
    take care of converting input values to the type of the corresponding
    column, and raise an exception if the conversion cannot be performed.
    '''

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

    @property
    def system(self):
        raise AttributeError('system')

    def keys(self):
        ''' sorted list of available properties '''
        return sorted(self.table.props)

    def __setitem__(self, prop, val):
        ''' update the value of prop with val '''
        p=self._ptr
        col=p.propIndex(prop)
        if col==BadId:
            raise KeyError, "No such property '%s'" % prop
        p.setProp(self._id, col, val)

    def __getitem__(self, prop):
        ''' get the value of prop '''
        p=self._ptr
        col=p.propIndex(prop)
        if col==BadId:
            raise KeyError, "No such property '%s'" % prop
        return p.getProp(self._id, col)

    def duplicate(self):
        ''' create a new entry in the parent parameter table with the
        same values as this one, returning it. '''
        return Param(self._ptr, self._ptr.duplicate(self._id))

class ParamTable(object):
    '''
    The `ParamTable` class is a 2d table, whose rows are indexed by ``id``
    and whose columns are properties; see the discussion of properties in
    the Overview.  A `ParamTable` is used by `TermTables` to hold the shared
    parameters for its `Terms`.
    '''

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

    def addParam(self, **kwds):
        ''' add and return a new Param().
        
        If keyword arguments are supplied, they will be assigned to the
        newly created Param before returning it.'''
        p=Param(self._ptr, self._ptr.addParam())
        for k,v in kwds.items():
            p[k]=v
        return p

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
    '''
    A `Term` is a handle for an entry in a `TermTable`.  

    The properties of a `Term` can be read and updated using a dictionary
    like interface.  Both "term properties" and properties from the
    `ParamTable` are accessed through the same interface.  To add or
    remove properties, use the provided methods in the `TermTable` or
    `ParamTable` instance.  If a `Term`'s ``param`` is shared by another
    `Term` in any other `TermTable`, Msys will take care of providing
    the `Term` with its own `Param` containing a copy of the original
    properties before applying the changes.  However, if you a modify a
    `Param` through its dictionary interface, you will affect all `Terms`
    that happen to share that `Param`::

        # fetch the stretch_harm table
        table = mol.table('stretch_harm')
        # update the properties of just the first Term
        table.term(0)['fc'] = 320
        # update the properties of all terms that use this param!
        table.term(0).param['fc'] = 320
    '''

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
    def system(self):
        ''' parent System of parent TermTable '''
        return TermTable(self._ptr).system

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
        if col==BadId:
            raise KeyError, "No such property '%s'" % prop
        id=self.paramid
        if id==BadId:
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
        if col==BadId:
            raise KeyError, "No such property '%s'" % prop
        id=self.paramid
        if id==BadId:
            raise RuntimeError, "No assigned param for '%s'" % repr(self)
        if p.refcount(id) > 1:
            id=p.duplicate(id)
            self.paramid=id
        p.setProp(id, col, val)

class TermTable(object):
    '''
    Each TermTable is intended to describe a specific type of interaction,
    e.g. stretch, angle, Lennard-Jones, constraint_hoh, etc.  A TermTable
    has an arity (given by the natoms property) which specifies how many
    atoms are involved in each interaction: one for nonbonded terms, two
    for stretch terms, etc.  Each interaction instance is described by
    a `Term`.  Each Term references the appropriate number of atoms,
    and exactly one Param, which lives in a ParamTable owned (or possible
    shared) by the TermTable.  

    The functional form described by a TermTable is not part of msys; all
    msys does is represent the forcefield parameters in a generic way.
    '''

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
    def props(self):
        ''' Table properties '''
        return PropertyMap(self._ptr.tableProps())

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
        ptr, ids = _convert_ids(atoms)
        if ptr!=self.system._ptr:
            raise ValueError, "atoms are from a different System"
        return [self.term(x) for x in self._ptr.findWithAll(ids)]

    def findWithAny(self, atoms):
        ''' return the terms that contain at least one of the given atoms '''
        if not atoms: return []
        ptr, ids = _convert_ids(atoms)
        if ptr!=self.system._ptr:
            raise ValueError, "atoms are from a different System"
        return [self.term(x) for x in self._ptr.findWithAny(ids)]

    def findWithOnly(self, atoms):
        ''' return the terms that contain only the given atoms '''
        if not atoms: return []
        ptr, ids = _convert_ids(atoms)
        if ptr!=self.system._ptr:
            raise ValueError, "atoms are from a different System"
        return [self.term(x) for x in self._ptr.findWithOnly(ids)]

    def findExact(self, atoms):
        ''' return the terms that contain precisely the given atoms in the 
        given order. '''
        if not atoms: return []
        ptr, ids = _convert_ids(atoms)
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
    '''

    The `System` class holds all structure and forcefield data
    for a single chemical system.  Create a new `System` using
    ``msys.CreateSystem()``, or from a file using ``msys.LoadDMS`` or
    ``msys.LoadMAE.``


    A `System` organizes the information in a DMS file into several
    different groups:

     * Tables - `TermTables` are grouped and accessed by name

     * cell - the unit cell vectors for the `System`, in the form of a 3x3
       NumPy array.

     * nonbonded_info - the NonbondedInfo object describing the type of
       nonbonded interactions.

     * provenance - a list of Provenance objects describing how the input
       file has been processed.

     * Auxiliary tables: Everything else in the DMS file that does not
       fit into one of the above categories finds its way into an auxiliary table.
       Notable denizens of this category include:

       - cmap tables

       - forcefield (annotation for parameters in the DMS file)

    '''

    __slots__ = ('_ptr', '_atoms')

    def __getinitargs__(self):
        ''' Pickle support (requires cPickle.HIGHEST_PROTOCOL) '''
        return self._ptr

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

    def hash(self):
        ''' hash of contents of this system.

        The hash is insensitive to provenance and System.name.
        However, adding atoms, bonds, terms, etc., then removing them may not
        leave the hash unchanged.
        '''
        return _msys.HashSystem(self._ptr)

    def asCapsule(self):
        ''' Return a capsule wrapper of the internal SystemPtr.

        The capsule holds a bare pointer and therefore must not outlive self.
        '''
        return self._ptr.asCapsule(self._ptr)

    @classmethod
    def fromCapsule(cls, cap):
        ''' Construct from a capsule wrapper of a SystemPtr.
        '''
        return cls(_msys.SystemPtr.fromCapsule(cap))

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

    def addChain(self, ct=BadId):
        ''' add and return a new Chain.
        If no ct is given, the chain will be added to the first ct,
        creating one if necessary. '''
        return Chain(self._ptr, self._ptr.addChain(ct))

    def addCt(self):
        ''' add and return a new Ct '''
        return Ct(self._ptr, self._ptr.addCt())

    def _remove(self, elems, klass):
        for a in elems:
            if a._ptr != self._ptr:
                raise ValueError("elements come from a different system")
            if not isinstance(a, klass):
                raise TypeError("expected list of %s; got %s" % (klass,type(a)))
            a.remove()

    def delAtoms(self, atoms):
        ''' remove the given Atoms from the System '''
        self._remove(atoms, Atom)

    def delBonds(self, bonds):
        ''' remove the given Bonds from the System '''
        self._remove(bonds, Bond)

    def delResidues(self, residues):
        ''' remove the given Residues from the System '''
        self._remove(residues, Residue)

    def delChains(self, chains):
        ''' remove the given Chains from the System '''
        self._remove(chains, Chain)

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

    def ct(self, id):
        ''' return the Ct with the specified id '''
        return Ct(self._ptr, id)

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
    @nonbonded_info.setter
    def nonbonded_info(self, nbinfo):
        self._ptr.nonbonded_info = nbinfo

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
    def ncts(self):
        ''' number of Cts '''
        return self._ptr.ctCount()

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
    def cts(self):
        ''' return list of all cts in the system '''
        ptr = self._ptr
        return [Ct(ptr, i) for i in ptr.cts()]

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
            self.cell[i][:] = cell[i]

    def getCell(self):
        ''' return copy of unit cell as 3x3 numpy array '''
        return self.cell.copy()

    @property
    def center(self):
        ''' return geometric center of positions of all atoms '''
        return numpy.mean(self.getPositions(), axis=0)

    def translate(self, xyz):
        ''' shift coordinates by given amount '''
        self._ptr.translate(*map(float, xyz))

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
    def topology(self):
        ''' list of bonded atoms for each atom in the System '''
        return self._ptr.topology()

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
        ''' Get the TermTable with the given name, raising ValueError if
        not present.
        '''
        ptr=self._ptr.table(name)
        if ptr is None:
            raise ValueError, "No such table '%s'" % name
        return TermTable(ptr)

    def getTable(self, name):
        ''' Return the TermTable with the given name, or None if not present.
        '''
        ptr=self._ptr.table(name)
        if ptr is None: return None
        return TermTable(ptr)

    def addTable(self, name, natoms, params = None):
        ''' add a table with the given name and number of atoms.
        If a table with the same name already exists, it is returned,
        otherwise the newly created table is returned.  If no ParamTable
        params is supplied, a new one is created.  '''
        if params is not None:
            params = params._ptr
        name = str(name)
        natoms = int(natoms)
        return TermTable(self._ptr.addTable(name, natoms, params))

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

    def atomsel(self, sel):
        ''' Create and return an atom selection object (Atomsel).
        Args:
            sel (object): str atom selection, or list of ids (possibly empty).

        Note:
            Even if ids are provided, the ids of the selection are in sorted order.
        '''
        if isinstance(sel, basestring):
            seltext = str(sel)
        elif not sel:
            seltext = 'none'
        else:
            seltext = 'index ' + ' '.join(map(str,sel))
        return Atomsel(self._ptr, seltext)

    def select(self, seltext):
        ''' return a list of Atoms satisfying the given VMD atom selection. '''
        ids=self._ptr.selectAsList(seltext)
        atms=self._update_atoms()
        return [atms[i] for i in ids]

    def selectIds(self, seltext, pos=None, box=None):
        ''' Return the ids of the Atoms satisfying the given VMD atom
        selection.  This can be considerably faster than calling select().

        if pos is supplied, it should be an Nx3 numpy array of positions,
        where N=self.natoms.

        If box is supplied, it should be a 3x3 numpy array of cell vectors,
        like System.cell.
        '''
        return self._ptr.selectAsList(seltext, pos, box)

    def selectArr(self, seltext):
        ''' Return the ids of the Atoms satisfying the given VMD atom
        selection as a numpy array of type uint32.  
        '''
        return self._ptr.selectAsArray(seltext)

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

    def selectCt(self, name=None):
        ''' Return a single Ct with the matching name, or raises an
        exception if no single such Ct is present '''
        cts = [c for c in self.cts if (name is None or c.name==name)]
        if not cts: return None
        if len(cts) is 1: return cts[0]
        raise ValueError, "Found %d cts with given name" % (len(cts))

    def append(self, system):
        ''' Appends atoms and forcefield from system to self.  Returns
        a list of of the new created atoms in self.  Systems must have
        identical nonbonded_info.vdw_funct.  Overwrites self.global_cell
        with system.global_cell only when self.global_cell is all zeros.
        '''
        p=self._ptr
        ids=p.append(system._ptr, BadId)
        atms=self._update_atoms()
        return [atms[i] for i in ids]

    def clone(self, sel=None, share_params=False):
        ''' Clone the System, returning a new System.  If selection is
        provided, it should be an atom selection string, a list of ids,
        or a list of Atoms.

        If share_params is True, then ParamTables will be shared between
        the old and new systems.  By default, copies of the ParamTables
        are made, but ParamTables shared _within_ the old system will
        also be shared in the new system.
        '''
        ptr = self._ptr
        if not sel:
            ids = ptr.atomsAsList()
        elif isinstance(sel, str):
            ids = ptr.selectAsList(sel)
        else:
            ids = list(sel)
            if isinstance(ids[0], Atom):
                ptr, ids = _convert_ids(sel)
                if ptr != self._ptr:
                    raise ValueError, "Atoms in sel are not from this System"
        flags = _msys.CloneOption.Default
        if share_params:
            flags = _msys.CloneOption.ShareParams

        return System(ptr.clone(ids, flags))

    def sorted(self):
        ''' Return a clone of the system with atoms reordered based on their 
        order of appearance in a depth-first traversal of the structure 
        hierarchy.
        '''
        return self.clone(self._ptr.orderedIds())

    def guessBonds(self, replace=True, reanalyze=True):
        ''' Guess bond connectivity based on an atomic-number based
        atom radius.  

        Replaces any existing bonds, unless replace=False is specified.

        Reanalyzes fragids and atom types unless reanalyze=False is specified.
        In that case, you MUST call updateFragids() manually before making
        any use of the fragment assignment (fragids will be out of date).
        '''
        if replace:
            self.delBonds(self.bonds)
        _msys.GuessBondConnectivity(self._ptr)
        if reanalyze:
            self.analyze()

    def analyze(self):
        ''' Assign atom and residue types.  This needs to be called
        manually only if you create a system from scratch, using 
        msys.CreateSystem(); in that case, analyze() should be called
        before performing any atom selections.
        '''
        _msys.Analyze(self._ptr)

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

    def findContactIds(self, cutoff,
            ids=None, other=None, pos=None):
        '''
        Find atoms not bonded to each other which are within cutoff of
        each other.  
        If ids is not None, consider only atoms with the given ids.  If
        other is not None, consider only atom pairs such that one is in ids
        and the other is in other.  If pos is not None, use pos as positions,
        which should be natoms x 3 regardless of the size of ids or other.
        pos may be supplied only when there are no deleted atoms in the
        structure.  

        Returns a list of (id 1, id 2, distance) tuples for each contact
        found.
        '''
        cutoff = float(cutoff)
        return self._ptr.findContactIds(cutoff, ids, other, pos)

class AnnotatedSystem(object):
    ''' System that has been annotated with additional chemical information

    The AnnotatedSystem class provides chemical annotation useful
    primarily for evaluating smarts patterns.  The system is expected to
    already have have chemical reasonable bond orders and formal charges,
    and to have no missing atoms (e.g. hydrogens).  If these criteria
    cannot be met, set allow_bad_charges=True in the constructor to bypass
    these checks; in that case the AnnotatedSystem can still be used to
    evaluate smarts patterns, but patterns making use of the electronic
    state of the system (e.g. aromaticity, hybridization, etc.) will
    not be correct (the system will appear to be entirely aliphatic).
    You may also use the AssignBondOrderAndFormalCharge function to
    assign reasonable bond orders and formal charges, assuming there
    are no missing atoms.

    The AnnotatedSystem defines a model for aromaticity.  First, the SSSR
    (smallest set of smallest rings) is determined.  Next, rings which
    share bonds are detected and grouped into ring systems.  Rings are
    initially marked as nonaromatic.  If the ring system taken as a whole
    is deemed to be aromatic, then all rings within it are aromatic as
    well; otherwise, individual rings are checked for aromaticity.  Rings
    are checked in this fashion until no new rings are found to be aromatic.
    
    A ring system is deemed to be aromatic if it satisfies Huckel's
    4N+2 rule for the number of electrons in the ring(s).  An internal
    double bond (i.e. a bond between two atoms in the ring) adds 2 to the 
    electron count.  An external double bond (a bond between a ring atom
    and an atom not in that ring) adds 1 to the electron count.  An
    external double bond between a carbon and a nonaromatic carbon makes
    the ring unconditionally nonaromtic.  An atom with a lone pair and
    no double bonds adds 2 to the electron count.
    '''

    def __init__(self, sys, allow_bad_charges=False):
        ''' Construct from System. AnnotatedSystem is not updated if System is
        subsequently modified. '''
        flags = _msys.AnnotatedSystemFlags.Default
        if allow_bad_charges:
            flags |= _msys.AnnotatedSystemFlags.AllowBadCharges
        self._ptr = _msys.AnnotatedSystem(sys._ptr, flags)
        self._name = sys.name

    def __repr__(self): return "<AnnotatedSystem '%s'>" % self._name

    @property
    def errors(self):
        ''' List of errors found during system analysis if 
        allow_bad_charges=True '''
        return self._ptr.errors()

    def aromatic(self, atom_or_bond):
        ''' Is atom or bond aromatic '''
        if type(atom_or_bond) == Atom:
            return self._ptr.atomAromatic(atom_or_bond.id)
        elif type(atom_or_bond) == Bond:
            return self._ptr.bondAromatic(atom_or_bond.id)
        else:
            raise TypeError, \
                    "atom_or_bond must be of type msys.Atom or msys.Bond"

    def hcount(self, atom):
        ''' Number of bonded hydrogens '''
        if type(atom) == Atom:
            return self._ptr.atomHcount(atom.id)
        else:
            raise TypeError, "atom must be of type msys.Atom"

    def degree(self, atom):
        ''' Number of (non-pseudo) bonds '''
        if type(atom) == Atom:
            return self._ptr.atomDegree(atom.id)
        else:
            raise TypeError, "atom must be of type msys.Atom"

    def valence(self, atom):
        ''' Sum of bond orders of all (non-pseudo) bonds '''
        if type(atom) == Atom:
            return self._ptr.atomValence(atom.id)
        else:
            raise TypeError, "atom must be of type msys.Atom"

    def loneelectrons(self, atom):
        ''' Number of lone electrons '''
        if type(atom) == Atom:
            return self._ptr.atomLoneElectrons(atom.id)
        else:
            raise TypeError, "atom must be of type msys.Atom"

    def hybridization(self, atom):
        ''' Atom hybridization -- 1=sp, 2=sp2, 3=sp3, 4=sp3d, etc.
        
        Equal to 0 for hydrogen and atoms with no bonds, otherwise
        max(1, a.degree() + (a.lone_electrons+1)/2 - 1). 
        '''
        if type(atom) == Atom:
            return self._ptr.atomHybridization(atom.id)
        else:
            raise TypeError, "atom must be of type msys.Atom"

    def ringbondcount(self, atom):
        ''' Number of ring bonds '''
        if type(atom) == Atom:
            return self._ptr.atomRingBonds(atom.id)
        else:
            raise TypeError, "atom must be of type msys.Atom"

class SmartsPattern(object):
    ''' A class representing a compiled SMARTS pattern 

    The Msys smarts implementation is similar to that of `Daylight smarts
    <http://www.daylight.com/dayhtml/doc/theory/theory.smarts.html>`, with
    support for arbitrarily nested recursive smarts.  A few features are 
    not currently supported; warnings will be generated when these constructs
    are used in a smarts pattern.
    
    * Directional bonds; e.g. ``\\`` and `/`; these are treated as single bonds
      (i.e. as a `-` character).
    
    * Chiral specification (``@``, ``@@``, etc); ignored.
    
    * Implicit hydrogen (``h``): treated as explicit ``H``.
    
    * Explicit degree (``D``): treated as bond count ``X``.
    
    * Isotopes: (``[12C]``): ignored.
    
    * Atom class (``[C:6]``): ignored.
    
    On the other hand, Msys does support hybridization using the ``^`` token,
    as in OpenBabel::
    
        [c^2]       select sp2 aromatic carbon

    '''
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

    def findMatches(self, annotated_system, atoms=None):
        ''' Return list of lists representing ids of matches of this pattern in
        this system, optionally requiring that the first atom match belongs to
        the given set of atoms. An AnnotatedSystem must be used here, which can
        be constructed from a System after calling
        AssignBondOrderAndFormalCharge. '''
        ptr = annotated_system._ptr
        if atoms is None:
            atoms = ptr.atoms()
        else:
            atoms = _convert_ids(atoms)[1]
        return self._pat.findMatches(ptr, atoms)

    def match(self, annotated_system):
        ''' Return True if a match is found anywhere; False otherwise.

        This is much faster than checking for an empty result from
        findMatches.
        '''
        return self._pat.match(annotated_system._ptr)

def CreateSystem():
    ''' Create a new, empty System '''
    return System(_msys.SystemPtr.create())

def _convert_ids(hs, klass=Atom):
    if not hs:
        return None, []
    ids = list()
    ptr = None
    for h in hs:
        if not isinstance(h, klass):
            raise TypeError("Expect list of %s; got %s" % (klass,type(h)))
        ids.append(h.id)
        ptr = ptr or h._ptr
        if ptr != h._ptr:
            raise ValueError("Inputs come from multiple systems")
    return ptr, ids

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
    of an MAE file, and the path argument will be ignored.

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

def Load(path, structure_only=False, without_tables=None):
    ''' Infer the file type of path and load the file.
    If without_tables is True or False, it specifies whether TermTables
    should be loaded along with the structure; by default it takes the
    value of structure_only.

    Returns a new System.
    '''
    structure_only = bool(structure_only)
    if without_tables is None:
        without_tables = structure_only
    ptr = _msys.Load(path, structure_only, without_tables)
    if not ptr:
        raise ValueError, "Could not guess file type of '%s'" % path
    return System(ptr)

class IndexedFileLoader(object):
    ''' Supports random access to multi-structure files
    '''

    def __init__(self, path, idx_path = None):
        ''' Open an indexed file loader, creating an index file if needed.
        Args:
            path (str): file path.  File type is inferred from the extension.
            idx_path (str): index file path.  Defaults to $path.idx.

        Note:
            You need write permission to the location of the index file.
        '''
        if idx_path is None:
            idx_path = ""
        else:
            idx_path = str(idx_path)
        self._ptr = _msys.IndexedFileLoader.create(path, idx_path)

    @property
    def path(self):
        ''' path to source file '''
        return self._ptr.path()

    def __len__(self):
        ''' number of entries '''
        return self._ptr.size()

    def __getitem__(self, index):
        ''' Get structure at given index

        Args:
            index (int): 0-based index

        Returns:
            mol (System): msys System
        '''
        return System(self._ptr.at(index))

def ConvertToRdkit(mol):
    ''' Construct an RDKit ROMol from the given System

    Args:
        mol (System): System

    Returns:
        rdkit.ROMol
    '''
    from rdkit import Chem
    rdmol = Chem.Mol()
    emol = Chem.EditableMol(rdmol)
    for atm in mol.atoms:
        ratm = Chem.Atom(atm.atomic_number)
        ratm.SetFormalCharge(atm.formal_charge)
        emol.AddAtom(ratm)
    for bnd in mol.bonds:
        emol.AddBond(bnd.first.id, bnd.second.id, Chem.BondType(bnd.order))
    rdmol = emol.GetMol()
    conf = Chem.Conformer(mol.natoms)
    for i, pos in enumerate(mol.getPositions()):
        conf.SetAtomPosition(i, pos)
    rdmol.AddConformer(conf)
    Chem.SanitizeMol(rdmol)
    Chem.AssignAtomChiralTagsFromStructure(rdmol)
    Chem.DetectBondStereoChemistry(rdmol, conf)
    Chem.AssignStereochemistry(rdmol)
    return rdmol

def LoadMany(path, structure_only=False, error_writer=sys.stderr):
    ''' Iterate over structures in a file, if the file type supports
    iteration.  

    for mol in LoadMany('input.mol2'): ...

    If there was an error reading a structure, LoadMany returns None for
    that iteration.  If error_writer is not None, it's write() method is
    invoked with the contents of the exception message as argument.  
    error_writer defaults to sys.stderr.
    '''
    it = _msys.LoadIterator.create(str(path), bool(structure_only))
    i=-1
    while True:
        i += 1
        try:
            mol = it.next()
        except Exception as e:
            if error_writer:
                error_writer.write("Error reading structure %d: %s\n" % (i,e.message))
            yield None
            continue
        if mol is None:
            break
        yield System(mol)

def Save(mol, path, append=False, structure_only=False):
    ''' Save the given system to path, using a file format guessed from the
    path name.  Not all formats support both append and structure_only options;
    see the corresponding SaveXXX functions.
    '''
    return _msys.Save(mol._ptr, str(path), _msys.Provenance.fromArgs(sys.argv),
            bool(append), bool(structure_only))

def FormatSDF(mol):
    ''' Return System in sdf format '''
    return _msys.FormatSDF(mol._ptr)

def ParseSDF(text):
    ''' Iterate over blocks in sdf format text.
    Accepts normal and gzipped text.
    '''
    it = _msys.ParseSDF(text)
    while True:
        mol = it.next()
        if mol is None:
            break
        yield System(mol)

def ReadPDBCoordinates(mol, path):
    ''' Read coordinates and box from the given pdb file into the given 
    System. 
    '''
    path=str(path)
    if not isinstance(mol, System): raise TypeError, "mol must be a System"
    _msys.ImportPDBCoordinates(mol._ptr, path)

def ReadCrdCoordinates(mol, path):
    ''' Read coordinates from the given Amber crd file into the given 
    System. 
    '''
    path=str(path)
    if not isinstance(mol, System): raise TypeError, "mol must be a System"
    _msys.ImportCrdCoordinates(mol._ptr, path)

def SaveDMS(system, path, structure_only=False, unbuffered=False):
    ''' Export the System to a DMS file at the given path. '''
    path = str(path)
    prov = _msys.Provenance.fromArgs(sys.argv)
    flags = _msys.DMSExportFlags.Default
    if structure_only:
        flags |= _msys.DMSExportFlags.StructureOnly
    if unbuffered:
        flags |= _msys.DMSExportFlags.Unbuffered
    if isinstance(system, System):
        _msys.ExportDMS(system._ptr, path, prov, flags)
    else:
        _msys.ExportDMSMany([x._ptr for x in system], path, prov)

def SerializeMAE(system, with_forcefield=True):
    ''' Return the MAE form of the System as a string. '''
    prov = _msys.Provenance.fromArgs(sys.argv)
    ff = bool(with_forcefield)
    flags = _msys.MaeExportFlags.Default
    if not with_forcefield:
        flags |= _msys.MaeExportFlags.StructureOnly
    return _msys.ExportMAEContents(system._ptr, prov, flags)

def SaveMAE(system, path, with_forcefield = True, append = False):
    ''' Export the System (or list of systems) to an MAE file at the 
    given path.  '''
    prov = _msys.Provenance.fromArgs(sys.argv)
    ff = bool(with_forcefield)
    flags = _msys.MaeExportFlags.Default
    if not with_forcefield:
        flags |= _msys.MaeExportFlags.StructureOnly
    if append:
        flags |= _msys.MaeExportFlags.Append
    
    if isinstance(system, System):
        _msys.ExportMAE(system._ptr, path, prov, flags)
    else:
        for mol in system:
            _msys.ExportMAE(mol._ptr, path, prov, flags)
            flags |= _msys.MaeExportFlags.Append

def SavePDB(system, path, append=False):
    ''' Export the System to a PDB file at the given path. '''
    flags = _msys.PDBExportFlags.Default
    if append:
        flags |= _msys.PDBExportFlags.Append
    _msys.ExportPDB(system._ptr, str(path), flags)

def SaveMol2(system, path, selection='none', append=False):
    ''' Export the System to a mol2 file at the given path.  You should
    probably call AssignBondOrderAndFormalCharge() before exporting
    the system.  '''
    ptr=system._ptr
    path=str(path)
    prov=_msys.Provenance.fromArgs(sys.argv)
    flags = _msys.Mol2ExportFlags.Default
    if append:
        flags |= _msys.Mol2ExportFlags.Append
    _msys.ExportMOL2(ptr,path,prov,system.selectIds(selection),flags)

def FromSmilesString(smiles, forbid_stereo=True):
    ''' Construct a System from a smiles string.

    Args:
        smiles (str): the smiles string
        forbid_stereo (bool): if True, raise exception if smiles has stereo

    EXPERIMENTAL.  In particular, stereo information in the smiles string
    is ignored.  Set forbid_stereo=False to permit stereo specifications
    to be silently ignored.  This flag may be removed at a later date once
    stereo support has been added.
    '''

    return System(_msys.FromSmilesString(smiles, forbid_stereo))

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
    ptr, ids = _convert_ids(atoms)
    rings = _msys.GetSSSR(ptr, ids, all_relevant)
    return [[Atom(ptr, id) for id in ring] for ring in rings]

def GetRingSystems(atoms):
    ''' Get ring systems for the given atoms '''
    ptr, _ids = _convert_ids(atoms)
    return _msys.RingSystems(ptr, _ids)

def AssignBondOrderAndFormalCharge(system_or_atoms, total_charge = None):
    """Assign bond orders and formal charges to a molecular system.

    Determines bond orders and formal charges by preferring neutral
    charges and placing negative charges with more electronegative
    atoms, under octet constraints and the total system charge
    constraint. Assigns the bond orders and formal charges to the system.
    Can assign to a subset of atoms of the system, provided these atoms
    form complete connected fragments.

    WARNING: calling this function on a chemically incomplete system,
    i.e. just protein backbone, may cause msys to hang indefinitely.

    Arguments:
    system_or_atoms: either a System or a list of Atoms
    total_charge: if not None, integral total charge
    """
    if isinstance(system_or_atoms, System):
        ptr = system_or_atoms._ptr
        if total_charge is None:
            return _msys.AssignBondOrderAndFormalCharge(ptr)
        ids = ptr.atomsAsList()
    else:
        ptr, ids = _convert_ids(system_or_atoms)

    if total_charge is None:
        _msys.AssignBondOrderAndFormalCharge(ptr, ids)
    else:
        _msys.AssignBondOrderAndFormalCharge(ptr, ids, int(total_charge))

class Graph(object):
    ''' Represents the chemical topology of a System

    Used mainly to implement graph isomorphism; see the match() method
    '''

    @classmethod
    def _from_boost(cls, _ptr, _sys):
        graph = cls(CreateSystem())
        graph._ptr = _ptr
        graph._sys = _sys
        return graph

    def __init__(self, system_or_atoms):
        if isinstance(system_or_atoms, System):
            ptr = system_or_atoms._ptr
            ids = ptr.atomsAsList()
        else:
            ptr, ids = _convert_ids(system_or_atoms)
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

def GuessHydrogenPositions(atoms):
    ''' Experimental '''
    ptr, ids = _convert_ids(atoms)
    _msys.GuessHydrogenPositions(ptr, ids)

def FindDistinctFragments(system):
    ''' Return fragids of representative fragments.  '''
    return _msys.FindDistinctFragments(system._ptr)

def ComputeTopologicalIds(system):
    ''' Compute and return the topological ids for the atoms or system '''
    ids = _msys.ComputeTopologicalIds(system._ptr)
    return [x for x in ids]

def CalcDistance(a, b):
    ''' Distance between atoms or positions a and b '''
    if not isinstance(a,numpy.ndarray): a = a.pos
    if not isinstance(b,numpy.ndarray): b = b.pos
    return _msys.calc_distance(a,b)

def CalcAngle(a, b, c):
    ''' Angle in radians of atoms or positions a-b-c. '''
    if not isinstance(a,numpy.ndarray): a = a.pos
    if not isinstance(b,numpy.ndarray): b = b.pos
    if not isinstance(c,numpy.ndarray): c = c.pos
    return _msys.calc_angle(a,b,c)

def CalcDihedral(a, b, c, d):
    ''' Dihedral angle in radians of atoms or positions a-b-c-d '''
    if not isinstance(a,numpy.ndarray): a = a.pos
    if not isinstance(b,numpy.ndarray): b = b.pos
    if not isinstance(c,numpy.ndarray): c = c.pos
    if not isinstance(d,numpy.ndarray): d = d.pos
    return _msys.calc_dihedral(a,b,c,d)

def ApplyDihedralGeometry(a, b, c, r, theta, phi):
    ''' Return the position of atom d with cd length r, bcd angle theta,
    and abcd dihedral phi, all in radians. 
    '''
    if not isinstance(a,numpy.ndarray): a = a.pos
    if not isinstance(b,numpy.ndarray): b = b.pos
    if not isinstance(c,numpy.ndarray): c = c.pos
    return _msys.apply_dihedral_geometry(a,b,c,float(r),float(theta),float(phi))

def CalcPlanarity(pos_or_atoms):
    ''' Planarity of positions or atoms '''
    if isinstance(pos_or_atoms, numpy.ndarray):
        return _msys.calc_planarity(pos_or_atoms)
    n = len(pos_or_atoms)
    pos = numpy.empty((n,3), 'd')
    for i in range(n):
        pos[i,:] = pos_or_atoms[i].pos
    return _msys.calc_planarity(pos)

def LineIntersectsTriangle(r, s, a, b, c):
    ''' True if line segment rs intersects triangle abc '''
    return _msys.line_intersects_tri(a,b,c,r,s)

class InChI(object):
    ''' InChI holds an the result of an inchi invocation for a structure '''

    def __init__(self, system, DoNotAddH = True, SNon = True, FixedH = True):
        opts = 0
        if DoNotAddH:
            opts |= _msys.InChI.Flags.DoNotAddH
        if SNon:
            opts |= _msys.InChI.Flags.SNon
        if FixedH:
            opts |= _msys.InChI.Flags.FixedH

        self._inchi = _msys.InChI.create(system._ptr, opts)

    def __str__(self):
        return self.string

    def __nonzero__(self):
        return self.ok

    @property
    def string(self):
        ''' Computed inchi string '''
        return self._inchi.string()

    @property
    def auxinfo(self):
        ''' Auxiliary info '''
        return self._inchi.auxinfo()

    @property
    def message(self):
        ''' Message returned by inchi during calculation '''
        return self._inchi.message()

    @property
    def ok(self):
        ''' Was an inchi computed successfully? '''
        return self._inchi.ok()

    @property
    def key(self):
        ''' inchi key for this object's string. '''
        return self._inchi.key()

def CloneSystem(atoms):
    ''' Call System.clone(atoms) using the System from the first atom.

    DEPRECATED.  Use System.clone directly instead.
    '''
    if not atoms: raise ValueError("empty atoms list")
    return atoms[0].system.clone(atoms)

class SpatialHash(object):
    ''' SpatialHash provides an interface for efficient spatial queries
    on particle positions. '''
    
    def __init__(self, pos, ids=None, box=None):
        ''' Construct from particle positions.  If ids are provided,
        they should be a numpy array of type uint32 and specify which
        rows of the Nx3 pos array are to be hashed.  
        If box is provided, it must be a 3x3 array of doubles, and the
        search will be performed using periodic boundary conditions.
        '''
        # FIXME: change the C++ SpatialHash class so that ids are optional.
        if ids is None:
            ids = numpy.arange(len(pos), dtype='uint32')
        self._hash = _msys.SpatialHash(pos, ids, box)

    def voxelize(self, radius):
        ''' Perform voxelization such that findWithin queries with 
        reuse_voxels=True at a radius equal to or less than the given 
        radius can be performed accurately.  For queries at radius 
        less than the voxelization, it may be worthwhile to revoxelize 
        at a smaller radius.  Note that, by default, findWithin calls
        voxelize with the query radius as argument, so it is not strictly
        necessary ever to use this method.
        '''
        self._hash.voxelize(float(radius))

    def findWithin(self, radius, pos, ids=None, reuse_voxels=False):
        ''' Find particles from pos which are within the given radius
        of some particle in the spatial hash (i.e. provided in the 
        SpatialHash constructor).  By default, voxelization is performed
        at the same resolution as the query radius, but this can be 
        overridden by calling voxelize() manually, then calling findWithin()
        with reuse_voxels=True.  pos is expected to be an Nx3 array of
        floats.  The ids parameter defaults to arange(len(pos)); supply
        an array of ids to limit the search to a subset of rows in pos.
        '''
        # FIXME: Make ids optional in the C++ interface
        if ids is None: ids = numpy.arange(len(pos), dtype='uint32')
        return self._hash.findWithin(radius, pos, ids, reuse_voxels)

    def findNearest(self, k, pos, ids=None):
        ''' Find at most k particles from pos with the smallest minimum 
        distance to some particle in the spatial hash.  If ids is not 
        provided, it defaults to arange(len(pos)).  
        '''
        if ids is None: ids = numpy.arange(len(pos), dtype='uint32')
        return self._hash.findNearest(k, pos, ids)

    def findContacts(self, radius, pos, ids=None, reuse_voxels=False):
        ''' Find pairs of particles within radius of each other.

        Args:
            radius (float): search radius
            pos (array[float]): positions
            ids (array[uint]): particle indices
            reuse_voxels (bool): assume voxelize(R>=radius) has already been called

        Returns:
           i, j, dists (tuple): Mx1 arrays of ids and distances.

        The first array corresponds to ids in the call to findContacts;
        the second column to the ids passed to the SpatialHash
        constructor.  

        IMPORTANT: pairs with the same id in the constructor and the call
        the findContacts are excluded from the output set.  Therefore,
        the positions passsed to findContacts should correspond to the
        same atom indices as the positions passed to the SpatialHash
        constructor.
        '''
        return self._hash.findContacts(radius, pos, ids, reuse_voxels)


HydrogenBond.__repr__ = lambda self: "<Hbond %s %s %s>" % (self.donor_id, self.acceptor_id, self.hydrogen_id)
HydrogenBond.donor = property(lambda x: x.donor_id, doc="Donor atom id")
HydrogenBond.acceptor = property(lambda x: x.acceptor_id, doc="Acceptor atom id")

class HydrogenBondFinder(object):
    ''' Find candidate hydrogen bonds.

    More hbonds will be found than are "realistic"; further filtering
    may be performed using the energy attribute of the returned hbonds.
    A reasonable filter seems to be around -1.0 (more negative is
    stronger); i.e. energies greater than that are more likely than not
    to be spurious.

    The HydrogenBond class can also be used directly to compute hydrogen
    bond geometry and energies by supplying donor, acceptor and
    hydrogen positions.
    '''
    def __init__(self, system, donors, acceptors, cutoff=3.5):
        '''
        Args:
            system (System): msys system
            donors: selection string, list of ids, or list of Atoms
            acceptors: selection string, list of ids, or list of Atoms
            cutoff (float): distance cutoff for donor and acceptor

        Note:
            If Atoms are provided, they must be members of system.
        '''
        if not donors: raise ValueError("Empty donors argument")
        if not acceptors: raise ValueError("Empty acceptors argument")
        if cutoff<0.1: raise ValueError("Unreasonably small cutoff %s" % cutoff)

        self.system = system
        if isinstance(donors, str):
            donors = system.selectIds(donors)
        elif isinstance(donors[0], Atom):
            donors = [a.id for a in donors]
        if isinstance(acceptors, str):
            acceptors = system.selectIds(acceptors)
        elif isinstance(acceptors[0], Atom):
            acceptors = [a.id for a in acceptors]
        self.cutoff = cutoff

        # find hydrogens for each donor
        dh = dict()
        for don in donors:
            atm = system.atom(don)
            hyd = [h.id for h in atm.bonded_atoms if h.atomic_number==1]
            if hyd:
                dh[don] = hyd
        self.donors = numpy.array(dh.keys(), dtype=numpy.uint32)
        self.acceptors = numpy.array(acceptors, dtype=numpy.uint32)
        self.hydrogens_for_donor = dh

    def find(self, pos=None):
        ''' Find hydrogen bonds for the given positions, defaulting to the
        current positions of the input system.'''
        if pos is None:
            pos = self.system.getPositions()
        elems = self.system.findContactIds(self.cutoff,
                                           self.donors,
                                           self.acceptors,
                                           pos)
        results = []
        hdict = self.hydrogens_for_donor
        for don, acc, dist in elems:
            hlist = hdict[don]
            if len(hlist) > 1:
                # use the hydrogen closest to the acceptor
                apos = pos[acc]
                hdist = [(CalcDistance(pos[h],apos), h) for h in hlist]
                hdist.sort()
                hyd = hdist[0][1]
            else:
                hyd = hlist[0]
            hbond = _msys.HydrogenBond(pos[don], pos[acc], pos[hyd])
            hbond.donor_id = don
            hbond.acceptor_id = acc
            hbond.hydrogen_id = hyd
            results.append(hbond)

        return results

