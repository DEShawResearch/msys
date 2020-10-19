import numpy
import os


class Atomsel(object):
    """ Supports alignment of molecular structures """

    __slots__ = ("_ptr", "_ids", "_seltext")
    """ atom selection object """

    def __init__(self, ptr, seltext):
        """ don't use directly - use System.atomsel() """
        self._ptr = ptr
        self._ids = ptr.selectAsArray(seltext)
        self._seltext = seltext

    def __len__(self):
        """ number of selected atoms """
        return len(self._ids)

    def __str__(self):
        return self._seltext

    def __repr__(self):
        return "<Atomsel '%s'>" % self._seltext

    @property
    def ids(self):
        """ ids of selected atoms in the parent system """
        return self._ids

    @property
    def system(self):
        """ parent system """
        from msys import System

        return System(self._ptr)

    def getPositions(self):
        return self._ptr.getPositions(self._ids)

    def _positions(self, other):
        if isinstance(other, Atomsel):
            opos = other.getPositions()
        elif isinstance(other, numpy.ndarray):
            opos = other
            if len(other) == self._ptr.atomCount():
                opos = other[self._ids]
            else:
                opos = other.copy()
        else:
            raise TypeError("Require either msys.Atomsel or numpy.ndarray")
        if len(self) != len(opos):
            raise ValueError(
                "Size mismatch: self (%d) != other (%d)" % (len(self), len(opos))
            )
        return opos

    def raw_alignment(self, other):
        """Compute alignment to other object.  Compute and return
        aligned rmsd, and rotational and translational transformations."""
        rpos = self.getPositions()
        opos = self._positions(other)
        if len(rpos) == 0:
            raise RuntimeError(
                "Empty atom selection '%s' - cannot compute alignment" % self
            )
        rcenter = rpos.mean(0)
        ocenter = opos.mean(0)
        rpos -= rcenter
        opos -= ocenter
        c = numpy.dot(rpos.transpose(), opos)
        v, s, w_tr = numpy.linalg.svd(c)
        is_reflection = (numpy.linalg.det(v) * numpy.linalg.det(w_tr)) < 0.0
        if is_reflection:
            s[-1] = -s[-1]
            v[:, -1] = -v[:, -1]
        rx = numpy.dot(v, w_tr).transpose()
        rpos *= rpos
        opos *= opos
        E0 = rpos.sum() + opos.sum()
        msd = (E0 - 2.0 * s.sum()) / len(rpos)
        msd = max([msd, 0.0])

        tx = rcenter - numpy.dot(ocenter, rx)
        return numpy.sqrt(msd), rx, tx

    def currentRMSD(self, other):
        """compute RMS distance to other object, which may be
        Atomsel or an array of positions.  In either it must be
        the case that len(other) equals len(self) or len(self.system)
        """
        dx = self.getPositions()
        dx -= self._positions(other)
        dx *= dx
        return numpy.sqrt(dx.sum() / len(dx))

    def alignedRMSD(self, other):
        """ Return the aligned rmsd to other. """
        return self.raw_alignment(other)[0]

    def alignCoordinates(self, other):
        """If other is an Atomsel instance, align the coordinates of
        other's System with self.  If other is a numpy array, align
        the array with self, using corresponding indices.

        In either case, return the aligned RMSD.
        """
        dx, rx, tx = self.raw_alignment(other)
        if isinstance(other, numpy.ndarray):
            pos = other
            pos[:] = numpy.dot(pos, rx)
            pos += tx
        elif isinstance(other, Atomsel):
            pos = other.system.getPositions()
            pos[:] = numpy.dot(pos, rx)
            pos += tx
            other.system.setPositions(pos)
        else:
            assert False

        return dx
