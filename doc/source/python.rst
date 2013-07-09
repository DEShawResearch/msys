````````````````
Python Scripting
````````````````

Most of the functionality in msys is exposed in its Python interface.
This section introduces the Python interface and explains how to use
it effectively.  We begin by introducing some concepts that pervade
the Python interface, then move to some examples.

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

Getting started
===============

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


Msys also lets you append chemical systems to an existing file, for certain
file formats.  The supported Save methods will have an 'append' option in
their function signatures.

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
  terms = stretch.terms
  params = stretch.params
  props = params.props # ['fc', 'r0']
  print "%d stretch terms, %d stretch params" % (terms.nterms, params.nparams)

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


  table = mol.addTable('funky_harm', 2)
  table.params.addProp('fk', float)
  table.params.addProp('r0', float)


If a table with a given name already exists in a `System`, ``addTable``
and ``addTableFromSchema`` will just return the existing table.



Files with multiple components
==============================

To examine every structure in a multi-component file without having to
load them all into memory at once, use **LoadMany**.  Unlike the **Load**
function, which always returns one System, **LoadMany** is a generator
which iterates over molecular structures in the input file::

  for mol in msys.LoadMany('input.mol2'):
     print mol.name

Not every file format supports LoadMany; in cases where it doesn't, LoadMany
will stop after a single iteration, yielding just one `System`.

If you use LoadMany to load a file, each `System` will have only one
`Ct`.  However, if you use Load to import an MAE or DMS file, and the
file contains multiple components, the new `System` will contain `Ct`
elements corresponding to those components::

  mol = msys.Load('small_vancomycin_complex.mae')
  for ct in mol.cts:
     print ct.name

  # prints:
  # vancomycin_diala_complex
  # SPC water box

The ct information wil be preserved when saving the System back to an MAE
or DMS file.  

You can create a multi-ct system from existing `Systems` using the
``append`` method::

  pro = msys.Load('protein.dms')
  pro.ct(0).name = 'protein'
  wat = msys.Load('water.dms')
  wat.ct(0).name = 'water'
  pro.append(wat)
  assert pro.ncts == 2     # assuming there was 1 ct in protein.dms and wat.dms
  assert pro.ct(1).name == 'water'
  msys.Save(pro, 'combined.dms')


