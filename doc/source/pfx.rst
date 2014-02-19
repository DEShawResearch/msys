
The Pfx module
===============================

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

TODO
----

 * Implement msys-style clone().  This might be needed for zendo, though
   it's not clear to me that zendo _really_ needs this.

 * Make the Pfx object serializable.  Again, might be needed for zendo,
   but maybe not.

 * Support triclinic cells, maybe.  We have no use for them internally.

 * Support weighted RMSD.


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

