
'''
This is the high-level Python interface for msys, intended for use
by chemists.
'''

import _msys

from _msys import GlobalCell, NonbondedInfo

class Handle(object):
    __slots__ = ('_ptr', '_id')

    def __init__(self, _ptr, _id):
        self._ptr = _ptr
        self._id = _id

    def __eq__(self, x): 
        try:
            return self._id==x._id and self._ptr==x._ptr
        except AttributeError:
            return False

    def __ne__(self, x): 
        try:
            return self._id!=x._id or  self._ptr!=x._ptr
        except AttributeError:
            return False

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
        return (self.x, self.y, self.z)
    @pos.setter
    def pos(self, xyz):
        self.x, self.y, self.z = xyz

    @property
    def vel(self):
        ''' velocity '''
        return (self.vx, self.vy, self.vz)
    @vel.setter
    def vel(self, xyz):
        self.vx, self.vy, self.vz = xyz

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
        'chargeB',
        'atomic_number', 'formal_charge',
        'moiety', 'alchemical',
        'name')

class Residue(Handle):
    __slots__ = ()

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
    def center(self):
        ''' return geometric center of positions of atoms in residue '''
        n=self.natoms
        if n==0: return [0,0,0]
        x=0.0
        y=0.0
        z=0.0
        for a in self.atoms:
            x += a.x
            y += a.y
            z += a.z
        return [x/n, y/n, z/n]

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

__add_properties(Chain, 'name')


class Param(object):
    __slots__=('_ptr', '_id')

    def __init__(self, _ptr, _id):
        self._ptr = _ptr
        self._id  = _id

    def __eq__(self, x): return self._id==x._id and self._ptr==x._ptr
    def __ne__(self, x): return self._id!=x._id or  self._ptr!=x._ptr
    def __repr__(self): return '<Param %d>' % self._id

    @property
    def id(self): 
        ''' id in parent table '''
        return self._id
    
    @property
    def table(self): 
        ''' parent ParamTable '''
        return ParamTable(self._ptr)

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

    def __eq__(self, x): return self._ptr==x._ptr
    def __ne__(self, x): return self._ptr!=x._ptr
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

    @property
    def shared(self):
        ''' is this ParamTable used by more than one TermTable? '''
        return self._ptr.shared()

class Term(object):
    __slots__=('_ptr', '_id', '_state')

    def __init__(self, _ptr, _id, _state='A'):
        self._ptr = _ptr
        self._id  = _id
        self._state = _state

    def __eq__(self, x): 
        return self._id==x._id and self._ptr==x._ptr and self._state==x._state
    def __ne__(self, x): 
        return self._id!=x._id or  self._ptr!=x._ptr or  self._state!=x._state

    def __repr__(self): return '<Term %d>' % self._id

    @property
    def paramid(self):
        if self._state is 'A':
            return self._ptr.param(self._id)
        return self._ptr.paramB(self._id)
    @paramid.setter
    def paramid(self, id):
        if self._state is 'A':
            self._ptr.setParam(self._id, id)
        else:
            self._ptr.setParamB(self._id, id)

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
    def paramA(self):
        ''' param for the A state of this Term '''
        return self.stateA.param
    @paramA.setter
    def paramA(self, val):
        self.stateA.param = val

    @property
    def paramB(self):
        ''' param for the B state of this Term '''
        return self.stateB.param
    @paramB.setter
    def paramB(self, val):
        self.stateB.param = val

    @property
    def alchemical(self):
        ''' Does this term contain an alchemical parameter set? '''
        return not _msys.bad(self._ptr.paramB(self._id))

    @property
    def stateA(self):
        ''' Returns the non-alchemical version of this Term '''
        if self._state=='A': return self
        return Term(self._ptr, self._id, 'A')

    @property
    def stateB(self):
        ''' Returns the alchemical (B) version of this Term '''
        if self._state=='B': return self
        return Term(self._ptr, self._id, 'B')

    @property
    def state(self):
        ''' The alchemical state this Term represents; either 'A' or 'B'. '''
        return self._state

    @property
    def atoms(self):
        ''' list of Atoms for this Term '''
        return [Atom(self._ptr.system(), i) for i in self._ptr.atoms(self._id)]

    @property
    def table(self): 
        ''' parent TermTable '''
        return TermTable(self._ptr)

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
        if ptr.paramRefs(id) > 1:
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
        self._ptr.system().removeTable(self._ptr)

    def coalesce(self):
        ''' Reassign param and paramB for each Term in this Table to a member
        of the distinct set of Params used by those Terms.
        '''
        self._ptr.coalesce()

    @property
    def name(self):
        ''' name of this table '''
        return self._ptr.name()
    @name.setter
    def name(self, newname):
        self._ptr.system().renameTable(self._ptr.name(), newname)

    def __repr__(self):
        return "<TermTable '%s'>" % self.name

    @property
    def params(self): 
        ''' The ParamTable for terms in this table. '''
        return ParamTable(self._ptr.params())

    @property
    def system(self): 
        ''' The System whose atoms are referenced by this table. '''
        return System(self._ptr.system())

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
        assert atom.system == self.system
        self._ptr.delTermsWithAtom(atom.id)

    @property
    def alchemical(self):
        ''' does this table contain alchemical terms? '''
        return self._ptr.alchemical()

    @property
    def terms(self):
        ''' returns a list of all the Terms in the table '''
        return [Term(self._ptr, i) for i in self._ptr.terms()]

    def term(self, id):
        ''' returns the Term in the table with the given id '''
        return Term(self._ptr, id)

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

class System(object):

    __slots__ = ('_ptr',)

    def __init__(self, _ptr):
        ''' Construct from SystemPtr.
        Do not invoke directly; use CreateSystem() instead.
        '''
        self._ptr = _ptr

    def __eq__(self, x): return self._ptr == x._ptr
    def __ne__(self, x): return self._ptr != x._ptr
    def __hash__(self): return self._ptr.__hash__()

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

    @property
    def positions(self):
        ''' Nx3 list of lists of positions of all atoms '''
        return self._ptr.getPositions()
    @positions.setter
    def positions(self, pos):
        self._ptr.setPositions(pos)

    @property
    def center(self):
        ''' return geometric center of positions of all atoms '''
        n=self.natoms
        if n==0: return [0,0,0]
        xyz=reduce(lambda a,b: (a[0]+b[0],a[1]+b[1],a[2]+b[2]), self.positions)
        return [a/n for a in xyz]

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

    def select(self, seltext):
        ''' return a list of Atoms satisfying the given VMD atom selection. '''
        p=self._ptr
        ids=p.select(seltext)
        return [Atom(p,i) for i in ids]

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

def file_endings(*args):
    ''' decorator: adds file_endings as list of case-insensitive regex patterns
    '''
    import re
    def _tmp(fcn):
        fcn.file_endings = [re.compile(".+\.%s$" % a, re.I) for a in args]
        return fcn
    return _tmp

@file_endings("dms", "dms.gz")
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


@file_endings("mae", "mae.gz", "maegz", "maeff", "maeff.gz", "cms", "cms.gz")
def LoadMAE(path=None, ignore_unrecognized = False, buffer=None):
    ''' load the MAE file at the given path and return a System containing it.
    Forcefield tables will be created that attempt to match as closely as
    possible the force terms in the MAE file; numerical differences are bound
    to exist.  If ignore_unrecognized is True, ignore unrecognized force
    tables.

    If the buffer argument is provided, it is expected to hold the contents
    of a DMS file, and the path argument will be ignored.

    If the contents of the file specified by path, or the contents of buffer,
    are recognized as being gzip-compressed, they will be decompressed on
    the fly. '''

    if buffer is None and path is None:
        raise ValueError, "Must provide either path or buffer"
    if buffer is not None and path is not None:
        raise ValueError, "Must provide either path or buffer"

    if path is not None:
        ptr = _msys.ImportMAE(path, ignore_unrecognized )
    else:
        ptr = _msys.ImportMAEFromBuffer( buffer, ignore_unrecognized )
    return System(ptr)

@file_endings("pdb")
def LoadPDB(path):
    ''' Load a PDB file at the given path and return a System.
    No bonds will be created, even if CONECT records are parent.
    '''
    path=str(path)
    ptr = _msys.ImportPDB(path)
    return System(ptr)

@file_endings("prmtop", "prm7")
def LoadPrmTop(path):
    ''' Load an Amber7 prmtop file at the given path and return a System.
    Coordinates and global cell information are not present in the file.
    '''
    path=str(path)
    ptr = _msys.ImportPrmTop(path)
    return System(ptr)

def GetLoader(path, defval=None):
    ''' Retrieve a procedure from the current module if it claims
    to be able to load files with the given path. If no such procedure is
    found, return defval '''
    for name, func in globals().items():
        for pat in getattr(func, "file_endings", ()):
            if pat.match(path):
                return func
    return defval

def Load(path):
    ''' Infer the file type of path using GetLoader and load the file.
    Returns a new System.
    '''
    func = GetLoader(path)
    if func is None:
        raise ValueError, "Could not guess file type of '%s'" % path
    return func(path)

def ReadCrdCoordinates(mol, path):
    ''' Read coordinates from the given Amber crd file into the given 
    System. 
    '''
    path=str(path)
    if not isinstance(mol, System): raise TypeError, "mol must be a System"
    _msys.ImportCrdCoordinates(mol._ptr, path)

def SaveDMS(system, path):
    ''' Export the System to a DMS file at the given path. '''
    import sys
    _msys.ExportDMS(system._ptr, path, _msys.Provenance.fromArgs(sys.argv))


def SaveMAE(system, path, with_forcefield = True ):
    ''' Export the System to an MAE file at the given path. '''
    import sys
    _msys.ExportMAE(system._ptr, str(path), bool(with_forcefield) )

def TableSchemas():
    ''' available schemas for System.addTableFromSchema '''
    return [s for s in _msys.TableSchemas()]

def NonbondedSchemas():
    ''' available nonbonded schemas for System.addNonbondedFromSchema '''
    return [s for s in _msys.NonbondedSchemas()]


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




