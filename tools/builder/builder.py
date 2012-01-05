
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

    def buildChain(self, chain):
        ''' build a particular chain of the given system using the previously
        loaded topology definitions '''
        assert isinstance(chain, msys.Chain)
        _builder.build(self._ptr, chain._ptr, chain.id)


def build(mol, defs):
    ''' build all chains of the given System using defs.  Return the
    new system with atoms put into natural order. '''
    assert isinstance(mol, msys.System)
    assert isinstance(defs, Defs)

    mol=mol.clone()
    for c in mol.chains: defs.buildChain(c)
    return mol.sorted()

