``````
User's Guide
``````

Introduction
============

Msys is a library for manipulating molecular structures and their
associated forcefields.  If you are already familiar with the Ent
library, you'll see that the two have much in common.

All molecular structure in Msys is held in an object called a `System`.
Within a `System`, individual particles, including physical atoms as
well as pseudo-particles used as interaction sites, are represented
as `Atoms`.  Bonds between `Atoms` are represented by a `Bond` 
structure.  `Atoms` are grouped into `Residues`, and `Residues`
are grouped into `Chains`.

A `System` also holds a set of `TermTables` representing the interactions
between `Atoms`.  A `TermTable` can be thought of as a particular kind
of interaction; for example, a fully parameterized system would likely
contain a ``stretch_harm`` `TermTable` to represent two-body covalent
bond forces.   Each `Term` in a `TermTable` refers to the same number
of atoms, though there can be any number of `Terms` in a given `TermTable`.

Typically, many of the interactions in a `TermTable` are parameterized
using identical parameters, especially when there are many identical
copies of the same molecule in the `System`.   For compactness, and also
for ease of forcefield parameterization, a `TermTable` holds a separate
table called a `ParamTable` which contains the interaction properties that
can be shared by many `Terms`.  Changes to an entry in a `ParamTable` will
affect the interaction strengths of every `Term` referencing that entry.
It is also possible for developers to construct multiple `TermTables`
that share the very same `ParamTable`, so that changes to a shared
`ParamTable` affect multiple `TermTables` or `Systems`.

In addition to this data representation, Msys also contains an atom
selection interface which implements most of the VMD atom selection
language.  It also comes with conversion tools and importers for
the Schrodinger MAE file format, and for the DESRES DMS file format.
The DMS file format is preferred due to its greater flexibility and
consistency.

Msys is written in C++, with bindings for Python.  The direct Python
bindings are generated using boost::python, and are intended to be
used only by developers, or those familiar with the C++ interface.
The rest of this guide describes the high level Python interface,
which is based on the Python bindings.  


Getting started
===============


