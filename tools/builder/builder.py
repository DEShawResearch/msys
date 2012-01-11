
import _builder
import msys

CHARMM='CHARMM'

class Defs(object):
    __slots__ = ('_ptr')

    def __init__(self):
        self._ptr = _builder.Defs()

    def load(self, path, format=CHARMM):
        ''' import topology definitions from the given file path. '''
        path=str(path)
        if format==CHARMM:
            self._ptr.import_charmm_topology(path)
        else:
            raise ValueError, "Unrecognized format '%s'" % format

    def buildChain(self, chain, pfirst="", plast=""):
        ''' build a particular chain of the given system using the previously
        loaded topology definitions '''
        assert isinstance(chain, msys.Chain)
        _builder.build(self._ptr, chain._ptr, chain.id, pfirst, plast)

    def patch(self, pres, *targets):
        ''' apply the patch residue with the given name to targets given
        by the remaining arguments. '''
        ids=msys._msys.IdList()
        mols=set()
        for t in targets:
            assert isinstance(t, msys.Residue)
            mols.add(t._ptr)
            ids.append(t.id)
            if len(mols)>1:
                raise ValueError, "all targets must be from the same System"
        if not mols:
            raise ValueError, "need at least one target"
        _builder.patch(self._ptr, pres, mols.pop(), ids)

    def build(self, mol):
        ''' build all chains of the given System using defs.  Return the
        new system with atoms put into natural order. '''
        assert isinstance(mol, msys.System)
    
        mol=mol.clone()
        for c in mol.chains: self.buildChain(c)
        return mol.sorted()

