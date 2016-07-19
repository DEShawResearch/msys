****************
Python Scripting
****************

Most of the functionality in msys is exposed in its Python interface.

Overview
========

This section introduces the Python interface and explains how to use
it effectively.  We begin by introducing some concepts that pervade
the Python interface, then move to some examples.

Msys ids
--------

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
---------------

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
---------------

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
---------------------------

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
------------------------------

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

The msys module
===============

.. automodule:: msys
    :members:
    :special-members:


Dealing with duplicate parameters
---------------------------------

After performing various modifications to a `TermTable`, you may find
that the associated `ParamTable` contains many entries whose values
are all identical.  The redundant parameters can be removed by first
"coalescing" the parameter assignments of each `Term` to a set of distinct
`Params`, then cloning the `System`.  When a `System` is cloned, only the
`Params` which are referenced by at least one `Term` in the `TermTable` are
copied to the new `System`:: 

  import msys
  mol=msys.CreateSystem()
  a1=mol.addAtom()
  a2=mol.addAtom()
  a3=mol.addAtom()
  table = mol.addTableFromSchema('stretch_harm')
  p1=table.params.addParam()
  p1['fc']=320
  p1['r0']=1.0
  t1=table.addTerm([a1,a2], p1)
  t2=table.addTerm([a1,a3], p1)
  
  # At this point we have two terms and one param.  Suppose we ignore the
  # fact that t1 and t2 share a Param, and we just update their properties
  # to the same value:
  
  t1['r0']=1.2
  t2['r0']=1.2
  
  # Now we have two Params, because when we updated t1, we created a second
  # Param that was unshared by t2.  When we updated t2, p1 was unshared, so
  # no duplicate was made.
  assert table.params.nparams==2
  
  # But we could get by with only a single Param.  Let's do that:
  mol.coalesceTables()
  
  # At this point t1 and t2 are sharing a Param, and the other one is unused:
  assert t1.param==t2.param
  assert table.params.nparams==2
  assert table.nterms==2
  
  # When we clone, the unused params are not copied to the new system.
  mol2=mol.clone()
  assert mol2.table('stretch_harm').params.nparams==1


.. autoclass:: msys.TermTable
   :members:


Sharing ParamTables
-------------------


Forcefield developers will (we hope!) appreciate the ability for Msys to
parameterize multiple `TermTables` from potentially different `Systems`
using a single `ParamTable` instance.  Normally, when a `System` is loaded
from an input file, or a `TermTable` is created using the scripting interface,
each `TermTable` refer to a `ParamTable` of its very own, and no other
`TermTable` can or will reference it.  However, at the time that a `TermTable`
is created, a `ParamTable` can be provided which will be used to hold
the `Param` entries for the `Terms` in the `TermTable`::


  # create two independent systems
  m1=msys.CreateSystem()
  m2=msys.CreateSystem()

  # add some atoms
  m1.addAtom()
  m2.addAtom()
  m2.addAtom()

  # create a free-standing ParamTable and add some Params
  params=msys.CreateParamTable()
  p1=params.addParam()
  p2=params.addParam()

  # create a table in system 1 which uses the free ParamTable
  table1=m1.addTable("table", 1, params)

  # no other TermTable is using the ParamTable
  assert not params.shared

  # create a table in system 2 which also uses the free ParamTable
  table2=m2.addTable("table", 1, params)

  # now the ParamTable is shared
  assert params.shared
  assert table1.params == table2.params

  # Add some terms to each table
  t1=table1.addTerm(m1.atoms, p2)
  t2=table2.addTerm(m2.atoms[1:], p2)

  assert t1.param == t2.param
  assert t2.param == p2

  # modifications to the the original table and its params are propagated
  # to each table
  params.addProp("fc", float)
  p1['fc']=32
  p2['fc']=42
  assert t1['fc']==42
  assert t2['fc']==42

  # p1 is shared by multiple TermTables, but within a TermTable, p1 is not
  # shared.  Modifications to t1['fc'] will affect t2!
  t1['fc'] = 52
  assert t2['fc'] == 52
  
Pfx
===

.. toctree::
   :maxdepth: 2

.. automodule:: msys.pfx
    :members:

What pfx does
-------------

**Pfx** can be configured to perform a number of tasks related to
postprocessing of trajectories of molecular systems.  There are four
main issues which **Pfx** is designed to deal with:

1. *Fixing bonds*.  If two atoms with a bond between them are found in
   different periodic images, one of them must be shifted by some integer
   linear combination of the unit cell vectors so that the distance
   between them is no greater than half a unit cell vector along each
   cell vector direction.  When there are multiple atoms bonded together,
   the bond fixing operation must be applied to the bonds composing this 
   connected component in topologically sorted order.

2. *Gluing components*.  Some molecular systems, such as multimeric
   ion channels, contain components which are not explicitly bonded
   to each other, but which do stay together during the simulation and
   should therefore be kept together during postprocessing.  For each
   set of glued components, **Pfx** finds the transformations which 
   minimize the square distance between the centers of each component.

3. *Centering and alignment*.  **Pfx** can either center a selected set
   of atoms on the origin, or align a selected set to a reference
   structure.  

4. *Wrapping components*.  Any of the preceeding operations could place
   the center of a connected set of atoms outside the unit cell
   centered at the origin.  **Pfx** shifts each connected component to
   bring it as close to the origin as possible, which maintaining any
   glued components.  Importantly, if an alignment has been performed
   in the previous step, then the rotational part of the alignment
   transformation must be applied to the unit cell before performing
   the wrapping.  Another subtlety is that when alignemnt has been
   performed, the wrapping should be performed about the center of the
   reference selection, not necessarily the origin.  Otherwise, if the
   reference structure is far from the origin, wrapping could undo the
   alignment.

The main work of **pfx** is does in the *apply* method.  The arguments
to *apply* are a coordinate set and, optionally, a periodic cell and/or
a set of velocities.  Here's what happens when you call *apply*, assuming
that both the periodic cell and the velocities have been provided:

#. Fix bonds.
#. Glue components.
#. Translate the entire system to bring the centered or aligned atoms to
   the origin.
#. Compute a rotational transformation which aligns the system to the
   centered reference structure.
#. Apply the rotation to the input positions, unit cel, and velocities.
#. Wrap connected and glued components.
#. Shift the entire system to the center of the reference structure.

Specifying topology
-------------------

A **Pfx** instance is constructed from a bond topology.  The bond topology
indicates both how many atoms are in the molecular system as well as
which atoms are bonded together.  **Pfx** analyzes this topology to
find the connected components comprising the system.  If the *fixbonds*
argument is True, then **Pfx** also computes and caches a topologically
sorted list of bonds from the topology, so that the bond fixing step
can be performed efficiently.


Specifying glue
---------------

The *glue* method of **Pfx** lets you specify atoms which should be kept
together even if there is no explicit bond between them.  Suppose the
ids of all the protein atoms in a multimeric protein are passed as the
argument to *glue*.  **Pfx** first finds the set of connected components
which overlap with the selection.  In the *glue components* step, the
centers of these components will be brought together.  Moreover, in the
*wrap components* step, all the protein atoms will be treated as a single
component for the purpose of wrapping.

Suppose now that only one residue from each monomer of a multicomponent
protein is included in a *glue* selection.  The same set of connected
components will be kept together as before, when the entire protein was
glued; however, the centers of the connected components will be computed
from just the glued residue in each monomer, rather than from all atoms
of each monomer.  The *wrap components* step will be unchanged.

The *glue* method can be called multiple times on the same **Pfx** instance.
It is perfectly valid for glue selections in different invocations to 
overlap.  


Performing both centering and alignment
---------------------------------------

When viewing trajectories, chemists often want to specify both "center"
and "fit" selections.  But what does this mean?  If you center on, say,
atoms 1-10, and align atoms 10-20, one operation will undo the other.
The only sensible approach seems to be to apply the "center" specification
to whatever is being used for the reference structure, and then use the
"fit" selection to align non-reference frames to that selection.


What about periodicfix?
-----------------------

The algorithms in this module are essentially the same as those in 
periodicfix.  So why a new module?  Here are some reasons:

 * Periodicfix isn't consistent about its use of float and double, and
   does a lot of interconversion.  This makes it slow.

 * Periodicfix doesn't have both single and double precision versions
   available from Python.  This one does.

 * Periodicfix makes you specify weights to align a subset of atoms. 
   Pfx doesn't use weights; or, if you like, the weights must be zero
   or one.  In practice it's been found that that's all anyone needs.
   Having to specify weights is cumbersome.  If someone really wants to
   have weight support in pfx we can add it some day.

 * Periodicfix has separate topology and fragment wrapper types, which
   make the Python interface more cumbersome.  Pfx has just one type.

 * Periodicfix has accreted additional functionality which has nothing
   to do with periodic images or alignment, including contact finding
   and a hydrogen bond energy function. 

 * The svd in periodicfix is greatly inferior to the one here.  Yes,
   it would be easy replace the one in periodicfix with this one.


Molfile
=======


.. automodule:: msys.molfile
    :members: Plugin, DtrReader, Frame, Atom, SeqFile, Grid
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

    .. method:: frame(n)

       returns the nth frame (0-based index).  Currently only the dtr
       plugin supports this method.

    .. method:: grid(n)

       return the nth grid.  For dx and ccp4 files.

Writer
------

.. autoclass:: msys.molfile.Writer
    :members:

    Writers are initialized with a path and either an array of Atoms or
    an atom count.  If the Writer supports structure writing, Atoms must
    be provided; if the Writer only writes frames, either one will do.

    .. method:: frame(f)

       If the writer supports frame writing, appends frame f to the end
       of the file.

    .. method:: grid(g)

       If the writer supports grid writing, writes Grid g to the file,
       where g is an instance of molfile.Grid, either returned from
       reader.grid(n) or created from scratch.

    .. method:: close()

       Invoked when the Writer goes out of scope, but it's not a bad
       idea to invoke it explicitly.


AnnotatedSystem
===============

.. autoclass:: msys.AnnotatedSystem
    :members:


SmartsPattern
=============

.. autoclass:: msys.SmartsPattern
    :members:


