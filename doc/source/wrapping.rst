
*****************
Periodic wrapping
*****************

Msys provides a high level python submodule named `wrap` for performing
periodic wrapping, and a low level module called `pfx`.

  
Pfx
===

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


Wrap
====

.. automodule:: msys.wrap
    :members:

