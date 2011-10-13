******
System
******

The `System` class holds all structure and forcefield data for a single
chemical system.  Create a new `System` using ``msys.CreateSystem()``, or 
from a file using ``msys.LoadDMS`` or ``msys.LoadMAE.`` 

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

NonbondedInfo
=============

.. autoclass:: msys.NonbondedInfo
   :members:


