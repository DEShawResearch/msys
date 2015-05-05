
import numpy as NP
from msys.pfx import Pfx

class Wrapper(object):
    def __init__(self, system, center=None, glue=[], unrotate=False):
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
        self.unrotate = unrotate

    def wrap(self):
        pos = self.mol.getPositions()
        box = self.mol.getCell()
        if self.unrotate:
            vecs = [r/NP.linalg.norm(r) for r in box]
            rot = NP.array(vecs).T
            pos = NP.dot(pos, rot)
            box = NP.dot(box, rot)
            self.mol.setCell(box)
        self.pfx.apply(pos,box)
        self.mol.setPositions(pos)

