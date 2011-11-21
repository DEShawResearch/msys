*********
TermTable
*********

Interactions between atoms in a `System` are described by the `TermTable`
class.  Each instance of the `TermTable` class describes a particular
kind of interaction; e.g., stretch terms, angle terms, dihedrals, etc.
Individual records within the `TermTable`, called `Terms`, describe the
interaction between particular sets of atoms.  Each `TermTable` holds one
`ParamTable`, which holds properties that can be shared by any number
of `Terms`.   A `TermTable` may also map additional properties, called
"term properties", to each `Term`.

Each `Term` in the TermTable uses the following attributes to represent
an interaction:

 * ``id``: the immutable id of the `Term`;

 * ``table``: the parent `TermTable`;

 * ``state``: the alchemical state; normally "A".  See below.

 * ``atoms``: a list of `Atoms`; the ``atoms`` in a `Term` cannot be changed;

 * ``param``: the `Param` in the `ParamTable` holding the parameters for
   the `Term`; multiple `Terms` may have the same ``param``.  ``param`` may
   be ``None``, indicating that the `Term` has not been assigned any 
   parameters.

The properties of a `Term` can be read and updated using a dictionary like
interface.  Both "term properties" and properties from the `ParamTable`
are accessed through the same interface.  To add or remove properties,
use the provided methods in the `TermTable` or `ParamTable` instance.
If a `Term`'s ``param`` is shared by another `Term` in the `TermTable`,
Msys will take care of providing the `Term` with its own `Param` containing
a copy of the original properties before applying the changes.  However,
if you a modify a `Param` through its dictionary interface, you will affect
all `Terms` that happen to share that `Param`::


  # fetch the stretch_harm table
  table = mol.table('stretch_harm')
  # update the properties of just the first Term
  table.term(0)['fc'] = 320
  # update the properties of all terms that use this param!
  table.term(0).param['fc'] = 320


We mentioned above that every `Term` has a ``state``, which is normally
"A".  In the "A" state, the ``param`` property corresponds to the
`Term` property ``paramA``.  However, every `Term` also has a ``paramB``
attribute, which is normally ``None``.  Setting ``paramB`` to a `Param`
from the `TermTable's` `ParamTable` makes the `Term` alchemical.  You can
access the properties of the alchemical state by fetching the "B" state
of the Term, using the ``stateB`` property.  A `Term` representing the
"B" state will read and write to its ``paramB`` member instead of ``paramA``, 
and the ``param`` property will correspond to ``paramB``.    


.. autoclass:: msys.TermTable
   :members:


Term
====

A `Term` is a handle for an entry in a `TermTable`.  

.. autoclass:: msys.Term
   :members:
   :special-members:
