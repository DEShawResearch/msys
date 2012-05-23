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
according to a combining rule.  In Msys, overrides to the parameters
in a `TermTable` are stored as a mapping from pairs of entries in the
``params`` to a entry in the ``override_params`` `ParamTable`.  Pairs
of `Params` are stored such that the ``id`` of the first `Param` is 
less than or equal to the ``id`` of the second `Param`; hence, there
are no redundant or conflicting overrides: if parameters *i* and *j*
have an override, then parameters *j* and *i* must be considered to
have the same override.


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


