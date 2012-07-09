
import numpy as NP
import periodicfix as PF

class Wrapper(object):
    def __init__(self, system, center=None, glue=[]):
        mol = system
        top = PF.Topology(mol.natoms)
        for b in mol.bonds:
            ai, aj = b.atoms
            top.add_bond(ai.id, aj.id)
        wrapper = PF.FragmentWrapper(top)
        if isinstance(glue, str):
            glue=[glue]
        for gsel in glue:
            for g in gsel.split(';'):
                sel = [a.id for a in mol.select(g)]
                wrapper.aggregate(sel)
        if center:
            center = [a.id for a in mol.select(center)]
            wrapper.aggregate(center)

        self.mol = mol
        self.top = top
        self.wrapper = wrapper
        self.center = center

    def wrap(self):
        pos = NP.zeros((self.mol.natoms, 3), 'f')
        box = NP.zeros((3,3), 'd')
        pos[:] = self.mol.positions
        for i in range(3):
            for j in range(3):
                box[i][j] = self.mol.cell[i][j]
        self.top.wrap_bonds(pos, box)
        self.wrapper.join(box, pos)
        if self.center:
            c=NP.zeros(3, 'd')
            for a in self.center:
                c += pos[a]
            c /= len(self.center) 
            pos -= c
        self.wrapper.wrap(box, pos, None)
        self.mol.positions = pos

