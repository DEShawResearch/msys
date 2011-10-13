*********
TermTable
*********

The `TermTable` class holds components of a molecular forcefield.  A given
`TermTable` instance belongs to a single System, and holds terms that
are all of the same type.  A `TermTable` instance be thought of as 
containing a number of `Terms`, which describe the interactions of
a single tuple of atoms.

Each `Term` in the TermTable is a handle contains three sorts of
attributes: (1) a fixed number of atoms, (2) parameters which are given
by the corresponding `Param` from its `ParamTable`, and (3) custom
"term properties" which are specific to the `Term` and not held by the
`ParamTable`.

The ``atoms`` in a Term cannot be changed once the `Term` is created; however,
its assigned `Param` can be set to any other `Param` from the `TermTable's`
own `ParamTable`, or to ``None``.

`Terms` can be added or removed from a `TermTable`; as with other handles,
the ``id`` of a `Term` never changes.  

.. autoclass:: msys.TermTable
   :members:


Term
====

A `Term` is a handle for an entry in a `TermTable`.  

.. autoclass:: msys.Term
   :members:
   :special-members:
