******
System
******

The `System` class holds all structure and forcefield data for a single
chemical system.  Create a new `System` using ``msys.CreateSystem()``, or 
from a file using ``msys.LoadDMS`` or ``msys.LoadMAE.`` 


A `System` organizes the information in a DMS file into several different
groups:

 * Tables - `TermTables` are grouped and accessed by name

 * cell - the GlobalCell object for the `System`.

 * nonbonded_info - the NonbondedInfo object describing the type of
   nonbonded interactions.

 * provenance - a list of Provenance objects describing how the input
   file has been processed. 

 * Auxiliary tables: Everything else in the DMS file that does not fit into
   one of the above categories finds its way into an auxiliary table.  
   Notable denizens of this category include:

   - cmap tables

   - forcefield (annotation for parameters in the DMS file)



.. autoclass:: msys.System
   :members:


Atom
====

.. autoclass:: msys.Atom
   :members:
   :special-members:

.. autoattribute:: msys.Atom.id
.. autoattribute:: msys.Atom.system

Bond
====

.. autoclass:: msys.Bond
   :members:
   :special-members:

.. autoattribute:: msys.Bond.id
.. autoattribute:: msys.Bond.system

Example::

  # create a list of bonded atoms for each atom in a System
  bondlists = [[b.id for b in a.bonded_atoms] for a in mol.atoms]
  # create a list of atom id pairs for each bond in a System
  bonds = [(b.first.id, b.second.id) for b in mol.bonds]

Residue
=======

.. autoclass:: msys.Residue
   :members:
   :special-members:

.. autoattribute:: msys.Residue.id
.. autoattribute:: msys.Residue.system

Chain
=====

.. autoclass:: msys.Chain
   :members:
   :special-members:

.. autoattribute:: msys.Chain.id
.. autoattribute:: msys.Chain.system


Ct
==

.. autoclass:: msys.Ct
   :members:
   :special-members:

.. autoattribute:: msys.Ct.id
.. autoattribute:: msys.Ct.system

GlobalCell
==========

The `GlobalCell` is a container for three vectors describing the periodic
nature of the `System`. 

.. autoclass:: msys.GlobalCell

   .. attribute:: A

      First shift vector

   .. attribute:: B

      Second shift vector

   .. attribute:: C

      Third shift vector

   .. method:: __getitem__(index)

      get corresponding shift vector.

    
Example::

  m=msys.CreateSystem()
  m.cell.A.x = 32           # A is now [32,0,0]
  m.cell[1][:] = [1,2,3]    # B is now [1,2,3]
  m.cell[2][1] = 5          # C is now [0,5,0]

