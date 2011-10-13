**********
ParamTable
**********

The `ParamTable` class is a 2d table, whose rows are indexed by ``id``
and whose columns are named and have a type of either `int`, `float`,
or `string`.  In Msys, these named, typed columns are referred to as
`properties`, or `props` for short.

A row in a `ParamTable` is queried and modified through a `Param` object,
described below.

As described in the Overview, a `ParamTable` is used by `TermTables`
to hold the shared parameters for its `Terms`.  However, properties are
used in several other places in Msys even when the underlying `ParamTable`
is not visible; for example, `Terms`, `Atoms`, and `Bonds` can all have
properties of their own, which are manipulated in much the same way
as properties in a `ParamTable`.  Familiarity with the `ParamTable`
interface will therefore make much of the rest of the Msys interface
easier to pick up.

.. autoclass:: msys.ParamTable
   :members:
   :special-members:

Param
=====

A `Param` instance is a reference to a row in a `ParamTable`.  Use the
``dict``-style interface to get and set values in the row.  Msys will
take care of converting input values to the type of the corresponding
column, and throw an error if the conversion cannot be performed.

.. autoclass:: msys.Param
   :members:
   :special-members:
