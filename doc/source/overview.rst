````````
Overview
````````

Introduction
============

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
However, as illustrated below, operations on an individual `Term` will affect
the interaction properties of just that `Term`; behind the scenes, Msys
takes care of creating a copy of a `Term`'s parameters as needed.

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

This section will guide you through some of the main features of Msys
in order to give you a feeling for how the Python interface works.

Once you have your Python environment by loading the appropriate
modules, fire up Python, import the msys module, and load a dms
or mae file::

  import msys

  # Load the entire contents of a DMS file
  dms=msys.LoadDMS('system.dms')

  # Import an MAE file, performing conversions on its forcefield data
  mae=msys.LoadMAE('system.mae')

You can also create a new `System` from scratch::

  # mol = msys.CreateSystem()

A `System` resides entirely in memory; changes to the `System` will not
persist until/unless you write it back out to a file::

  # Save the system as a DMS file
  msys.SaveDMS(dms, 'output.dms')

  # Export to MAE file
  msys.SaveMAE(dms, 'output.mae')


The full set of `Atoms`, `Bonds`, `Residues`, `Chains`, and `TermTables`
are available by fetching them from the system.   You can also fetch
the bonds involving a particular atom, the atoms in a residue, or the bonds
in a chain in a similar way::

  # get the number of atoms, and the total charge
  atoms = dms.atoms
  natoms = len(atoms)
  total_charge = sum(a.charge for a in atoms)

  # find the atoms participating in double bonds
  for chn in mol.chains:
    for res in chn.residues:
      for atm in res.atoms:
        for bnd in atm.bonds:
          if bnd.order == 2:
            print "atom %d in chain %s residue %s:%d has a double bond" % (
              atm.id, chn.name, res.name, res.num)

  # iterate over tables, print atoms per term and number of terms
  for t in mol.tables:
    print "table %s: %d atoms, %d terms" % (t.name, t.natoms, t.nterms)

  # fetch the stretch_harm table.  Throws an exception if no such table
  stretch = mol.table('stretch_harm')

Atom selections let you fetch a list of atoms using the VMD atom selection
language.  The ``select`` method returns a list of `Atoms`, which is
just a subset of the list that would be returned by the ``atoms`` property::

  # fetch the backbone atoms.  Note that bb is just a Python list
  bb = mol.select('backbone')


Once you have the atoms, if you actually want to work with
the residues or chains, it's easy to do::

  # get the set of distinct residues in the backbone
  bb_residues = set(a.residue for a in bb)

Note that the atoms returned by ``select`` refer back to the original
system.  Msys also provides the means to create a new `System` independent
of the original, using either the ``CloneSystem`` function or the 
``clone`` method of System.  When you clone a subset of a `System`, the 
`Terms` in the forcefield whose atoms are completely contained in the 
selection will be copied to the new `System`::

  # get the list of protein atoms
  pro_atoms = mol.select('protein')

  # construct a new system containing just the protein
  protein = msys.CloneSystem(pro_atoms)

  # Atoms in the cloned system have the same attributes as the originals,
  # but modifications to one do not affect the other
  assert pro_atoms[0].charge == protein.atoms[0].charge
  pro_atoms[0].charge += 3
  assert pro_atoms[0].charge != protein.atoms[0].charge

The ``clone`` method of `System` is a more concise way of selecting a
set of atoms, then immediately creating a new `System` from it::

  # create a new System with all the hydrogens removed
  hless = mol.clone('not hydrogen')

  # create a copy of the original
  dup = mol.clone()

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

You can change the properties of a selected `Term` using a 
dictionary-like interface::

  # Change the force constant of the first stretch term to 42
  stretch.terms[0]['fc] = 42


Adding new forcefield terms
===========================

Msys provides an interface for adding a `TermTable` corresponding
to a "standard" forcefield term and configuring that table with
its category and its the expected set of properties::


  # Get the available set of TermTable schemas:
  schemas = msys.TableSchemas()

  # For bonded, constraint, virtual, and polar terms, as well as 
  the exclusion table:
  table = mol.addTableFromSchema('posre_harm')  # position restraints

  # Get the available set of nonbonded schemas
  nb_schemas = msys.NonbondedSchemas()

  # For a nonbonded table:
  nb = mol.addNonbondedFromSchema('vdw_12_6')


The ``addNonbondedFromSchema`` also takes care of configuring the
``nonbonded_info`` properties of the `System`; see the section on
nonbonded parameters for more details.

If you have a new table type that hasn't made it into Msys' canonical
set, you can simply use ``addTable`` and configure the table yourself::


  table = mol.addTable('funky_harm')
  table.params.addProp('fk', float)
  table.params.addProp('r0', float)


If a table with a given name already exists in a `System`, ``addTable``
and ``addTableFromSchema`` will just return the existing table.


Msys properties
===============

Many objects in Msys (in particular, `Atoms`, `Bonds`, `Terms`, and
`Params`) can have typed attributes given to all members of the set
to which they belong.  In Msys, these attributes are referred to as
`properties`, or `props` for short, and have a type of either `int`,
`float`, or `str` (string).  The available property names and their
types can be queried in the appropriate parent object, using the
``props``, ``atom_props``, etc. properties of the parent.
The value of the property for a given element can be read and modified
using a dictionary-like interface on the element itself::

  mol = msys.LoadDMS('input.dms')
  # find all distinct values of the 'grp_energy' atom property, if it exists
  grp_energy_vals = set()
  if 'grp_energy' in mol.atom_props:
    for atm in mol.atoms:
      grp_energy_vals.add( atm['grp_energy'] )

  # add a new property 'foo' of type 'float'
  mol.addAtomProp('foo', float)
  # Set the value of foo to the z coordinate of the atom
  for a in mol.atoms: a['foo'] = a.pos[2]

When you add a property to a set of elements, the initial value will be 0
for `int` and `float` types, and the empty string for `str` types.  If a
property with the same name and type already exists, no action is taken.
An exception is thrown if you try to add a property with the same name 
but different type from an existing property.


Msys ids
========

In Msys, instances of the `Atom`, `Bond`, `Residue`, and `Chain` classes
are all `Handles`, in the sense that they refer to a piece of data held
by the parent `System`.  All Msys handles have an immutable ``id``
property that uniquely identifies them within their parent `System`.
Objects that hold references to other objects do so through the ``id``
of that object.  Two handles of the same type will compare equal to each
other if and only if they belong the same `System` and possess the same
``id``.

When you load a system from a file, or create one from scratch, these
``ids`` will be numbered consecutively, starting at zero.  Deleting
`Atoms`, `Bonds`, etc. from the `System` can introduce gaps in the set of
``ids``, but, once created, the ``id`` of an object never changes.

When Msys writes a DMS file, the primary keys of the particles will
be contiguous starting at 0, and will appear in the order in which the
particles appear in the `System`, even if the ``ids`` of the atoms in the
`System` are noncontiguous due to deletions.  When Msys loads a DMS file,
if the primary keys happen to be noncontiguous, Msys will still create a
`System` with the usual contiguous ids.

Notes on atom selections
========================

Msys implements essentially all of the atom selection language of VMD.
Differences between Msys and VMD's implementations (other than as yet
undiscovered bugs in Msys!) fall into the following categories:

* Element matching: In Msys, the atom selections "carbon", "hydrogen",
  "oxygen", etc. are based on the atomic number of the atoms.  In 
  VMD, you maybe shocked and surprised to learn, these atom selections
  are computed using a regular expression based on the atom name::
  
    vmd > atomselect macro oxygen
    name "O.*"
    
    vmd > atomselect macro hydrogen
    name "[0-9]?H.*"
    
    vmd > atomselect macro nitrogen
    name "N.*"
    
    vmd > atomselect macro carbon
    name "C.*" and not ion
  

  It was felt that, rather than slavishly follow VMD in this respect, Msys
  should try to get the correct answer.  Do you really want your "nitrogen"
  atom selection to include sodium (NA)?

* Field size: DMS and MAE files can hold chain, segment, and residue names
  of arbitrary length.  In Msys, these values are used as-is.  In VMD,
  the values are truncated; in particular, chain will be truncated to
  a single character in VMD, but not by Msys.

* Data representation: Msys has no concept of secondary structure, so the
  "sheet", helix", etc. atom selection keywords are not implemented in 
  msys.
  
* Floating-point roundoff: There may occasionally be differences in the
  results of distance based atom selections simply due the fact that Msys
  stores positions as doubles, while VMD stores them as floats.  

Notes on loading chemical systems
=================================

As mentioned earlier, Msys groups all `Atoms` into `Residues`, and all
`Residues` into `Chains`.  This hierarchy is, unfortunately, rarely made
explicit in the chemical system files in wide use, so Msys must infer the
grouping based on the values of certain particle attributes.

Msys uses the ``resname``, ``resid``, and ``chain`` names in the particle
table to groups atoms into `Residues`.  Any particle with the same value
for these three attributes will wind up in the same `Residue`, regardless of
the order of the particles or their bond structure.  In addition, starting
with Msys version 1.0.14, the ``segid`` particle attribute is also used
to disambiguate residues.  

Msys assigns particles to chains using only the ``chain`` attribute.  This
attribute can be of any length, although some legacy file formats restrict
it to only one character.  The ``segid`` attribute is not used for assigning
chain; if present in the file, it simply appears as an additional particle
attribute, and, in versions of Msys later than 1.0.14, will also be used to
disambiguate residues.

Whitespace in atom and residue names
------------------------------------

The PDB file format specifies that atom and residue names should be
aligned to particular columns within a 4-column region.  Unfortunately,
some have taken this alignment requirement to mean that an atom's
name actually includes the surrounding whitespace!  When Msys loads
a chemical system, atom and residue names are stripped of leading and
trailing whitespace before they are inserted into the structure.  


