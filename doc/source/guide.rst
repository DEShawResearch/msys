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

Once you have your Python environment by loading the appropriate
modules, fire up Python, import the msys module, and load a dms
or mae file::


  import msys
  dms=msys.LoadDMS('system.dms')
  mae=msys.LoadMAE('system.mae')

The full set of `Atoms`, `Bonds`, `Residues`, `Chains`, and `TermTables`
are available by fetching them from the system.   You can also fetch
the bonds involving a particular atom, the atoms in a residue, or the bonds
in a chain in a similar way::

  # get the number of atoms, and the total charge
  atoms = dms.atoms
  natoms = len(atoms)
  total_charge = sum(a.charge for a in atoms)

  # iterate over chains, then residues, then atoms:
  for chn in mol.chains:
    for res in chn.residues:
      for atm in res.atoms:
        for bnd in atm.bonds:
          a1, a2 = bnd.atoms
          bond_order = bnd.order

  # iterate over tables, print atoms per term and number of terms
  for t in mol.tables:
    print "table %s: %d atoms, %d terms" % (t.name, t.natoms, t.nterms)

  # fetch the stretch_harm table
  stretch = mol.table('stretch_harm')

Atom selections let you fetch a list of atoms using the VMD atom selection
language.  The ``atomselect`` method returns a list of `Atoms`, which is
just a subset of the list that would be returned by the ``atoms`` property::

  # fetch the backbone atoms.  Note that bb is just a Python list
  bb = mol.atomselect('backbone')


Once you have the atoms, if you actually want to work with
the residues or chains, it's easy to do::

  # get the set of distinct residues in the backbone
  bb_residues = set(a.residue for a in bb)

Note that the atoms returned by ``atomselect`` refer back to the original
system.  Msys also provides the means to create a new `System` independent
of the original, using the ``CloneSystem`` function::

  # get the list of protein atoms
  pro_atoms = mol.atomselect('protein')

  # construct a new system containing just the protein
  protein = msys.CloneSystem(pro_atoms)

  # Atoms in the cloned system have the same attributes as the originals,
  # but modifications to one do not affect the other
  assert pro_atoms[0].charge == protein.atoms[0].charge
  pro_atoms[0].charge += 3
  assert pro_atoms[0].charge != protein.atoms[0].charge

You can append the structure and associated forcefield from one `System`
onto another using System's ``append`` method::

  # duplicate the protein by appending to itself
  protein.append(protein)

  # load a water system and append it to the protein system.  Just as for
  # CloneSystem, after appending water to protein, modifications to water
  # will not affect any atoms in protein.
  water = msy.LoadDMS('water.dms')
  protein.append(water)

Terms in a system's forcefield can be accessed and modified by going 
through the corresponding `TermTable`::

  stretch = protein.table('stretch_harm')
  props = stretch.props # ['fc', 'r0']
  terms = stretch.terms
  params = stretch.param
  print "%d stretch terms, %d stretch params" % (len(terms), len(params))

Note that, because parameters can be and often are shared, modifications
to the properties of one `Term` can affect the properties of many others.
In the previous code snippet, it is likely that there will be many fewer
params than terms, indicating that some terms share params.  



