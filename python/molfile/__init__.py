"""
Structure and coordinate file manipulation library.

Reading a structure file::

    reader = molfile.mae.read('/path/to/foo.mae')

Iterating through the frames in a file::

    for frame in molfile.dtr.read('/path/to/foo.dtr').frames():
        function( frame.pos, frame.vel, frame.time, frame.box )

Random access to frames (only dtr files support this currently)::

    f27 = molfile.dtr.read('/path/to/foo.dtr').frame(27) # 0-based index

Write a trajectory to a frameset (dtr)::

    f = msys.molfile.Frame(natoms)
    w = msys.molfile.dtr.write('output.dtr', natoms=natoms)
    for i, xyz in enumerate(xyzs):
        f.pos[:] = xyz
        f.time = i
        w.frame(f)
    w.close()

Convert an mae file to a pdb file::

    input=molfile.mae.read('foo.mae')
    output=molfile.pdb.write('foo.pdb', atoms=input.atoms)
    output.frame(input.frames().next())
    output.close()

Write every 10th frame in a dtr to a trr::

    input=molfile.dtr.read('big.dtr')
    output=molfile.trr.write('out.trr, natoms=input.natoms)
    for i in range(0,input.nframes, 10):
        output.frame( input.frame(i) )
    output.close()


Write a frame with a specified set of gids::

    f = molfile.Frame(natoms, with_gids=True
    f.gid[:] = my_gids
    f.pos[:] = my_positions
    w.frame(f)

Read the raw fields from a frameset (dtr)::

    dtr = molfile.DtrReader('input.dtr')    # also works for stk
    for i in range(dtr.nframes):
        f = dtr.frame(i)
        keyvals = dict()
        frame = dtr.frame(i, keyvals=keyvals)
        ## use data in keyvals


Write raw fields to a frameset (dtr)::

    dtr = molfile.DtrWriter('output.dtr', natoms=natoms)
    keyvals = dict( s = "a string",
                    f = positions.flatten(),    # must be 1d arrays
                    i = numpy.array([1,2,3]),
                    )
    dtr.append( time = my_time, keyvals = keyvals )

"""

from ._molfile import *
from . import findframe
import os, sys
import bisect
import time
import numpy
import os

extensiondict = dict()


def register_plugin(plugin):
    """ put plugin in the global namespace, and add to extensiondict """
    globals()[plugin.name] = plugin
    d = extensiondict
    for ext in plugin.filename_extensions.split(","):
        ext = ext.strip()
        if ext:
            d.setdefault(ext, []).append(plugin)


_molfile.register_all(register_plugin)


def guess_filetype(filename, default=None):
    """ return plugin name based on filename, or default if none found. """
    dot = filename.rfind(".")
    if dot > 0:
        key = filename[dot + 1 :].lower()
        val = extensiondict.get(key)
        if val:
            return val[-1]
    return default


def list_fields(stk_or_tr):
    """list available fields in an stk or frameset (etr or dtr)"""
    trj = DtrReader(stk_or_tr)
    kv = dict()
    trj.frame(0, keyvals=kv)
    return sorted(kv.keys())


class FrameIter(object):
    def __init__(self, reader):
        if reader.nframes >= 0:
            self.reader = reader
            self.curframe = 0
            self.nframes = self.reader.nframes
        else:
            # number of frames is not known, so we have to use the
            # Reader.next() method which does not rewind.  Reload the
            # reader so that we can keep a private copy of our
            # current frame.
            self.reader = reader.reopen()
            self.curframe = -1

    def __iter__(self):
        return self

    def __next__(self):
        if self.curframe < 0:
            # iterate the reader itself
            f = self.reader.next()
            if not f:
                raise StopIteration
            return f

        # using our own iterator state
        if self.curframe >= self.nframes:
            raise StopIteration
        f = self.reader.frame(self.curframe)
        self.curframe += 1
        return f

    def next(self):
        return self.__next__()

    def skip(self, count=1):
        r = self.reader
        if count < 0:
            raise ValueError("skip count must be nonnegative")
        while count:
            r.skip()
            count -= 1


def reader_frames(self):
    return FrameIter(self)


_molfile.Reader.frames = reader_frames
del reader_frames


class Grid(object):
    def __init__(self, data, name="", axis=None, origin=None):
        """construct a new Grid object.
        data - 3d array of data; a copy is made.
        name - name for the grid.  default empty string.
        axis - 3x3 array with grid axes in the rows.  Default diag(1,1,1)
        origin - 0,0,0 corner of grid.  Defaults to [0,0,0]
        """
        self._data = numpy.array(data, "f")
        dims = self._data.shape
        if len(dims) != 3:
            raise ValueError("bad dims: want 3, got %d" % (len(dims)))

        self._name = name

        self._axis = numpy.eye(3, dtype="f")
        if axis is not None:
            axis = numpy.array(axis, "f")
            if axis.shape != (3, 3):
                raise ValueError("bad axis: want 3x3, got %s" % (axis.shape,))
            self._axis[:] = axis

        self._origin = numpy.zeros(3, "f")
        if origin is not None:
            origin = numpy.array(origin, "f")
            if origin.shape != (3,):
                raise ValueError("bad origin: want 3 elements, got %d" % (len(origin)))
            self._origin[:] = origin

    def __repr__(self):
        x, y, z = self._data.shape
        return "<%s: (%d, %d, %d)>" % (self.name, x, y, z)

    @property
    def name(self):
        return self._name

    @property
    def data(self):
        return self._data

    @property
    def axis(self):
        return self._axis

    @property
    def origin(self):
        return self._origin


def _grid_from_reader(reader, n):
    meta = reader.grid_meta(n)
    data = numpy.empty(meta["dims"], "f")
    reader.grid_data(n, data)
    name = meta["name"]
    x = meta["xaxis"]
    y = meta["yaxis"]
    z = meta["zaxis"]
    origin = meta["origin"]
    return Grid(data, name=name, axis=[x, y, z], origin=origin)


def _grid_to_writer(writer, grid):
    d = {
        "name": grid.name,
        "origin": grid.origin.tolist(),
        "xaxis": grid.axis[0].tolist(),
        "yaxis": grid.axis[1].tolist(),
        "zaxis": grid.axis[2].tolist(),
        "dims": grid.data.shape,
    }
    writer._grid(d, grid.data)
    return writer


_molfile.Reader.grid = _grid_from_reader
_molfile.Writer.grid = _grid_to_writer


class StkFile(object):
    """ Generalized stk file: handles any molfile format that provides times"""

    name = "ls"
    filename_extensions = "ls"

    @classmethod
    def read(cls, path, filetype=None):
        return cls.Reader(path, filetype)

    class Reader(object):
        def __init__(self, path, filetype):
            # path is location of stkfile
            # Extract the dirname
            # Is path a symlink? If so, figure out where the real file is.
            dirname = os.path.dirname(os.path.realpath(path))
            with open(path) as fp:
                lines = (x.strip() for x in fp)
                paths = [x for x in lines if x]
                paths = [
                    os.path.join(dirname, x) if not os.path.isabs(x) else x
                    for x in paths
                ]
            if filetype is None:
                klasses = (guess_filetype(path, SeqFile) for path in paths)
            else:
                klass = globals()[filetype]
                klasses = (klass for path in paths)
            readers = [klass.read(r) for klass, r in zip(klasses, paths)]

            # remove trailing empty readers
            while readers and readers[-1].nframes == 0:
                del readers[-1]
            # Store the list of sizes to allow get_prop to work.
            self._sizes = [r.nframes for r in readers]
            if readers:
                first = readers[-1].times[0]
                i = len(readers) - 1
                while i > 0:
                    i -= 1
                    cur = readers[i]
                    n = cur.nframes
                    while n > 0 and cur.times[n - 1] >= first:
                        n -= 1
                    self._sizes[i] = n
                    if n > 0:
                        first = min(first, cur.times[0])

            seqmap = []
            offset = 0
            for sz in self._sizes:
                offset += sz
                seqmap.append(offset)

            self.readers = readers
            self.seqmap = seqmap
            self.times = numpy.concatenate(
                [r.times[:sz] for r, sz in zip(self.readers, self._sizes)]
            )

        @property
        def natoms(self):
            return 0 if not self.readers else self.readers[0].natoms

        @property
        def nframes(self):
            return 0 if not self.readers else self.seqmap[-1]

        def frames(self):
            return (self.frame(i) for i in range(self.nframes))

        def frame(self, n):
            index = bisect.bisect_right(self.seqmap, n)
            if index > 0:
                n -= self.seqmap[index - 1]
            return self.readers[index].frame(n)

        def get_prop(self, prop):
            """
            Use the same technique as used for times to generate a properly ordered set of properties across all readers
            """
            return numpy.concatenate(
                [r._get_prop(prop)[:sz] for r, sz in zip(self.readers, self._sizes)]
            )

        def at_time_near(self, time):
            return self.frame(findframe.at_time_near(self.times, time))


register_plugin(StkFile)

with_pandas = True


class SeqFile(object):
    """Read csv-like files with column names in the first row"""

    filename_extensions = "seq"
    name = "seq"

    @classmethod
    def read(cls, path):
        """ Open an eneseq file for reading """
        return cls.Reader(path)

    class Reader(object):
        @staticmethod
        def _parse_header(line):
            # mapping from eneseq column name to molfile.Frame attribute
            ene2frame = dict(
                time="time",
                E="total_energy",
                E_p="potential_energy",
                E_k="kinetic_energy",
                E_x="extended_energy",
                P="pressure",
                V="volume",
                T="temperature",
            )
            elems = [x for x in line[1:].split() if ":" in x]
            if not elems:
                # handle metadynamics output.
                attrs = [x for x in line[1:].split()]
            else:
                attrs = (x.split(":")[1] for x in elems)
            return [ene2frame.get(x, x) for x in attrs]

        def __init__(self, path):
            """
            Only use a direct read to determine the header for the file
            Otherwise, use numpy.loadtxt
            """
            props = None
            with open(path) as fp:
                for line in fp:
                    cols = line.strip().split()
                    if len(cols) > 1 and cols[0] == "#" and cols[1].endswith("time"):
                        props = self._parse_header(line)
                        break

            self.rows = None
            if with_pandas:
                import pandas

                try:
                    self.rows = pandas.read_csv(
                        path, delim_whitespace=True, header=None, comment="#"
                    ).values
                except ValueError:
                    pass
            if self.rows is None:
                self.rows = numpy.loadtxt(path, ndmin=2)
            self.times = self.rows[:, 0] if len(self.rows) else numpy.array([])
            self.props = props

        @property
        def natoms(self):
            return 0

        @property
        def nframes(self):
            return len(self.rows)

        def _get_prop(self, prop):
            """
            Returns a MUTABLE reference to a property.
            Used to eliminate a copy when concatenating multiple seq files
            """
            try:
                return self.rows[:, self.props.index(prop)]
            except ValueError:
                raise ValueError("Cannot find property {:s}".format(prop))

        def get_prop(self, prop):
            # Column accessor
            return numpy.copy(self._get_prop(prop))

        def frame(self, n):
            f = Frame(natoms=0)
            for key, val in zip(self.props, self.rows[n]):
                setattr(f, key, val)
            return f

        def frames(self):
            return (self.frame(i) for i in range(self.nframes))

        def at_time_near(self, time):
            return self.frame(findframe.at_time_near(self.times, time))


register_plugin(SeqFile)
