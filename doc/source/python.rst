****************
Python Reference
****************

The msys module
===============

.. automodule:: msys
    :members:
    :special-members:
    :inherited-members:
    :exclude-members: __lt__, __gt__, __eq__, __ne__, __hash__, __getitem__, __new__, __repr__, __setitem__, __contains__, __delitem__

Molfile
=======


.. automodule:: msys.molfile
    :members: Plugin, DtrReader, Frame, Atom, SeqFile, Grid
    :inherited-members:
    :undoc-members:

Reader
------

.. autoclass:: msys.molfile.Reader
    :members:

    A Reader is a handle to an open file.  Use the atoms member to fetch the
    atomic structure from the file, assuming it exists.  To access frames,
    there are two methods.

    .. method:: frames()

       returns a FrameIter object for iteration over frames.  FrameIter
       has two methods: the usual next() method which returns a Frame,
       and skip(n=1), which advances the iterator by n frames without
       (necessarily) reading anything.  FrameIter is a very poor iterator:
       once a frame has been read or skipped, it can't be loaded again;
       you have use a brand new Reader.

    .. method:: grid(n)

       return the nth grid.  For dx and ccp4 files.

Writer
------

.. autoclass:: msys.molfile.Writer
    :members:

    Writers are initialized with a path and either an array of Atoms or
    an atom count.  If the Writer supports structure writing, Atoms must
    be provided; if the Writer only writes frames, either one will do.

    .. method:: grid(g)

       If the writer supports grid writing, writes Grid g to the file,
       where g is an instance of molfile.Grid, either returned from
       reader.grid(n) or created from scratch.


