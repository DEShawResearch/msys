**********
ParamTable
**********

The `ParamTable` class is a 2d table, whose rows are indexed by ``id``
and whose columns are properties; see the discussion of properties in
the Overview.  A `ParamTable` is used by `TermTables` to hold the shared
parameters for its `Terms`.

.. autoclass:: msys.ParamTable
   :members:
   :special-members:

Param
=====

A `Param` instance is a reference to a row in a `ParamTable`.  Use the
``dict``-style interface to get and set values in the row.  Msys will
take care of converting input values to the type of the corresponding
column, and raise an exception if the conversion cannot be performed.

.. autoclass:: msys.Param
   :members:
   :special-members:
