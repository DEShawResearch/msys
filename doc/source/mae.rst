
MAE file support
================

As of this writing (January 2013), the MAE file format still sees 
common use within DESRES and beyond.  The following notes are intended
to guide the expectations of users of MAE files when interacting
with msys and msys-derived tools such as viparr.  In particular,
we take note of MAE file contents that cannot be represented
within Msys at all, due to specific aspects of the design of Msys
that conflict with MAE files, and vice versa: not all Msys systems
can be readily serialized to an MAE file.

Multiple ct blocks
------------------

When there is more than one "f_m_ct" block in a MAE file, it can mean
either that the file contains multiple structures, or that the individual
blocks are separated for bookkeeping purposes but should be considered
to be a single physical system.  The Desmond and Anton software takes
the latter point of view, and it is therefore necessary to merge or
otherwise reconcile redundant information in the ct blocks when treating
them as a single system.

The LoadMany function in Python (LoadIterator in C++) loads each ct block
separately into a new System.

The LoadMAE function in Python (ImportMAE in C++) loads all ct blocks in
an MAE file and merges them into a single System.  The following 
requirements must be satisfied on the input blocks:

 * The global cell parameters (r_chorus_box_{abc}{xyz}) of the _last_
   ct block takes precedence.

 * The nonbonded functional form in ffio_vdwtypes, and the combining
   rule, must be the same in all blocks.  If the rule is left blank,
   rules in prior blocks takes precedence.

Note that the mapping of atoms to ct block is preserved.



Ct-level attributes
-------------------

When using LoadMany, the "name" attribute of each System will be read from
the "m_title" field of the corresponding ct block.  LoadMAE will set the
System name to the path name of the MAE file.

Msys reads and writes a ct-level array called "msys_provenance" in order
to track the history of its files.

Msys reads all the key-value pairs at the top level of the f_m_ct block, 
and the attributes in the m_depend subblock if it exists.  These attributes
are stored in the Ct structure of the Systems.

Pair terms
----------

Mae files offer several ways of representing pair terms in the ffio_pair
table.  The most common way is to store a scale factor for the LJ and
electrostatic components of the pair, though it is also possible to
store the pair cofficients explicitly.  Msys only stores pair terms with
their explicit interaction coefficients.  When reading MAE files with
scaled pair terms, the coefficients are calculated using the current
VDW parameters and particle charges.  When writing MAE files, Msys
always writes the explicit cofficients.


Pseudoparticles
---------------

Within a ct block, in addition to the physical atoms in the m_atom block,
there may be pseudoparticles in the ffio_pseudos block.  Msys has only
a single list of particles in a system, and thus must make a choice as
to the relative order of the atoms and pseudos.  The choice made is the
one also used by Desmond: within a ct block, all the atoms come first,
followed by all the pseudos.  In particular, the order of entries within
the ffio_sites block makes no difference in the ordering of the particles.

This means that if an Msys system contains atoms interleaved with pseudos,
as is often convenient, the system cannot be serialized to an MAE file.


VDW overrides (NBFix)
---------------------

NBFix terms are fully supported by Msys, which reads ffio_vdwtypes_combined
directly into the corresponding nonbonded Overrides table.

Multiple nonbonded tables
-------------------------

Although Msys permits multiple nonbonded tables to be defined, there is
currently no facility to represent these tables in MAE files.

Supported MAE tables
--------------------

Msys understands the following forcefield tables and functional forms
(ffio_funct field):

 * ffio_angles: ffio_funct = "ub" or "harm".  Urey-bradley "ub" terms
   are mapped to the "stretch_harm" table; "harm" terms are mapped
   to the "angle_harm" table.

 * ffio_bonds: ffio_funct = "harm"; all terms map to the "stretch_harm"
   table.

 * ffio_cmap[1-6]: ffio_ai, ffio_aj, and ffio_c1 fields map to 
   an auxiliary table with corresponding properties "phi", "psi", and
   "energy".

 * ffio_constraints: maps to "constraint_xyz" where "xyz" is the
   corresponding ffio_funct value.

 * ffio_dihedrals: ffio_funct = "proper_trig", "improper_trig",
   "opls_proper", or "opls_improper" map to the "dihedral_trig" table;
   ffio_funct = "improper_harm" maps to the "improper_harm" table.

 * ffio_exclusions: maps to the "exclusion" table.

 * ffio_{posre,angle,improper,stretch}_fbhw.

 * ffio_pairs: Supported vdw functional forms are "lj12_6_sig_epsilon", 
   "polynomial_cij", and "exp_6x".  Supported combining rules are
   "geometric", "arithmetic/geometric", and "lb/geometric".  Supported
   ffio_funct values are "coulomb", "coulomb_scale", "coulomb_qij",
   "lj", "lj_scale", and "lj_12_6_sig_epsilon".  Pair terms map to
   the "pair_12_6_es" table for "lj12_6_sig_epsilon" and "polynomial_cij"
   systems, and to the "pair_exp_6_es" table for "exp_6x" systems.

 * ffio_restraints: ffio_funct = "harm" maps to the "posre_harm" table.

 * ffio_torsion_torsion: maps to "torsiontorsion_cmap" table.

 * ffio_virtuals: Terms are mapped to tables named "virtual_xyz" where
   "xyz" is the ffio_funct value for the term.



