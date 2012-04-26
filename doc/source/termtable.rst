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

 * ``atoms``: a list of `Atoms`; the ``atoms`` in a `Term` cannot be changed;

 * ``param``: the `Param` in the `ParamTable` holding the parameters for
   the `Term`; multiple `Terms` may have the same ``param``.  ``param`` may
   be ``None``, indicating that the `Term` has not been assigned any 
   parameters.

The properties of a `Term` can be read and updated using a dictionary like
interface.  Both "term properties" and properties from the `ParamTable`
are accessed through the same interface.  To add or remove properties,
use the provided methods in the `TermTable` or `ParamTable` instance.
If a `Term`'s ``param`` is shared by another `Term` in any other `TermTable`,
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



Dealing with duplicate parameters
=================================

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
===================


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
  


Term
====

A `Term` is a handle for an entry in a `TermTable`.  

.. autoclass:: msys.Term
   :members:
   :special-members:

