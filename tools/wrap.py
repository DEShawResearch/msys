"""
Tools for periodic wraping
"""

import numpy as NP
from msys.pfx import Pfx
import msys


class Wrapper(object):
    """
    Example:

        wrap = msys.wrap.Wrapper(system, center="protein")
        wrap.wrap()
    """

    def __init__(
        self, system: msys.System, center: str = None, glue=[], unrotate=False
    ):
        """Construct a Wrapper for a system"""
        mol = system
        pfx = Pfx(mol.topology, fixbonds=True)
        if isinstance(glue, str):
            glue = [glue]
        for gsel in glue:
            for g in gsel.split(";"):
                pfx.glue(mol.selectIds(g))
        if center:
            pfx.align(mol.selectIds(center))

        self.mol = mol
        self.pfx = pfx
        self.unrotate = unrotate

    def wrap(self):
        """perform periodic wrapping using current coordinates"""
        pos = self.mol.getPositions()
        box = self.mol.getCell()
        if self.unrotate:
            vecs = [r / NP.linalg.norm(r) for r in box]
            rot = NP.array(vecs).T
            pos = NP.dot(pos, rot)
            box = NP.dot(box, rot)
            self.mol.setCell(box)
        self.pfx.apply(pos, box)
        self.mol.setPositions(pos)
