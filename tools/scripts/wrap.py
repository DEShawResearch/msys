
import numpy as NP
from msys.pfx import Pfx

class Wrapper(object):
    def __init__(self, system, center=None, glue=[]):
        mol = system
        pfx = Pfx(mol.topology, fixbonds=True)
        if isinstance(glue, str):
            glue=[glue]
        for gsel in glue:
            for g in gsel.split(';'):
                pfx.glue(mol.selectIds(g))
        if center:
            pfx.align(mol.selectIds(center))

        self.mol = mol
        self.pfx = pfx

    def wrap(self):
        pos = self.mol.getPositions()
        box = self.mol.getCell()
        self.pfx.apply(pos,box)
        self.mol.setPositions(pos)

