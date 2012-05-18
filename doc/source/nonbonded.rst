********************
Nonbonded parameters
********************


Nonbonded parameters for particles in Msys are handled as follows.
A `System` may have at most one `TermTable` whose name is "nonbonded" and
whose category is also "nonbonded".  The `ParamTable` for the nonbonded
table, not surprisingly, holds the nonbonded parameters for all the atoms.
Atoms are assigned a nonbonded type by creating `Terms` in the nonbonded
`TermTable`.   There should be exactly one `Term` for each `Atom`, and 
each `Atom` should be represented by exactly one `Term`.  

The `System` class has a method called ``addNonbondedFromSchema`` which is
a shortcut for creating a nonbonded table of a particular type.  The argument
to ``addNonbondedFromSchema`` will be the ``vdw_funct`` that appears in
the ``nonbonded_info`` field of the `System`.  The following ``vdw_funct``
values are currently supported:

 * ``vdw_12_6`` : parameters ``sigma``, ``epsilon``

 * ``vdw_exp_6`` : parameters ``alpha``, ``epsilon``, ``rmin``

 * ``vdw_exp_6s`` : parameters ``sigma``, ``epsilon``, ``lne``

When a `System` is created by loading a DMS or MAE file, a nonbonded
table will be created if and only if the input file contains nonbonded
information.  When saving a `System` to a DMS file, Msys checks that
there is at most one nonbonded table, and if one exists, ensures that
every `Atom` is found in exactly one `Term`.  

Overriding nonbonded interactions
=================================

Nonbonded interactions between particles are usually calculated by looking
up the nonbonded parameters (e.g., charge, sigma, epsilon) of the two
interacting particles, performing some sort of combining operation on
those parameters (e.g., geometric mean of the charge, arithmetic mean
of the sigma), then using those values in the functional form of the
interaction.

The DMS and MAE file formats allow one to specify nonbonded types whose
combined values are to be taken from a table, rather than computed
according to a combining rule.  In Msys, it is not possible to store
this type of table directly, because each `Term` in a `TermTable`
references only one `Param`, whereas an interaction override involves two.
Storing the combined tables as an auxiliary table isn't acceptable,
either, because any change to the particle nonbonded assignment would
silently break the intended override.

Instead, when Msys loads a "nonbonded_combined_param" table from a DMS
file, it creates a two-particle `TermTable` called "nonbonded_combined"
which holds the interaction override for each affect pair of particles.
The category of this table is "OVERRIDE", which lets Msys recognize
the table on DMS export and convert it back to its original form.

The 'msys.vdw.Combine' function implements the operation of overriding
the interaction between two sets of particles while leaving the other
interactions alone.  The 'dms-override-vdw' script wraps this operation
in a command line tool.

Alchemical nonbonded interactions
=================================

DMS files use a table called *alchemical_particle* to indicate which
particles have alchemical nonbonded states, and the parameters for
those states.  Msys represents the information in that table with
a `TermTable` called *alchemical_nonbonded*.   This table will share
a `ParamTable` with the regular *nonbonded* table, but will contain
`Terms` only for the alchemical particles.  The parameter for each
`Term` in *alchemical_nonbonded* will correspond to the B state of
the term's particle.  Additional per-particle information, such
as *chargeB*, *chargeC*, or *moiety*, will appear as term properties
for the particles.


NonbondedInfo
=============

.. autoclass:: msys.NonbondedInfo
   :members:


