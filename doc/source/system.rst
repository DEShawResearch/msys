******
System
******

The `System` class holds all structure and forcefield data for a single
chemical system.  Create a new `System` using ``msys.CreateSystem()``, or 
from a file using ``msys.LoadDMS`` or ``msys.LoadMAE.`` 


A `System` organizes the information in a DMS file into several different
groups:

 * Tables - `TermTables` are grouped and accessed by name

 * cell - the unit cell vectors for the `System`, in the form of a 3x3
   NumPy array.

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

