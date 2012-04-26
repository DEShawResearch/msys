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
according to a combining rule.  The way this is done is through an
auxiliary table named "nonbonded_combined_param".  This table should
have the same properties as the nonbonded table, and in addition
the integer properties "param1" and "param2", which refer to the id
of the nonbonded paramter being overridden.  

For example, suppose one wishes to zero out the vdw interaction between
atoms with ids 100 and 200.  Here is a script that will create a
new nonbonded combined table and override the interaction between just
those two particles::

  # load the system and fetch the nonbonded table
  mol=msys.LoadDMS('input.dms')
  nb=mol.table('nonbonded')

  # create the nonbonded combined table
  combined = msys.CreateParamTable()
  combined.addProp('param1', int)
  combined.addProp('param2', int)
  combined.addProp('sigma', float)
  combined.addProp('epsilon', float)

  # add the nonbonded combined table to the system
  mol.addAuxTable('nonbonded_combined_param', combined)

  # find the nonbonded parameters for the two atoms, and duplicate them.
  a1=mol.atom(100)
  a2=mol.atom(200)
  p1=None
  p2=None
  for t in nb.terms:
    if t.atoms[0]==a1:
      p1=t.param.duplicate()
      t.param=p1
    elif t.atoms[0]==a2:
      p2=t.param.duplicate()
      t.param=p2
  assert p1 is not None and p2 is not None

  # add a row to nonbonded combined for the these two parameters
  row = combined.addParam()
  row['param1'] = p1.id
  row['param2'] = p2.id
  row['sigma'] = 0.0
  row['epsilon'] = 1.0

  # paranoia: make the nonbonded combined table symmetric
  row = combined.addParam()
  row['param1'] = p2.id
  row['param2'] = p1.id # swapped p1 and p2
  row['sigma'] = 0.0
  row['epsilon'] = 1.0

  # save the result
  msys.SaveDMS(mol, 'combined.dms')

We duplicated the parameters for particles 100 and 200 because those
entries in the parameter table may be in use by other particles, and
we did not want to affect the interactions of any particles except
that pair.  In the new system, particles 100 and 200 will experience no
vdw interaction with each other, though they will experience the same
electrostatic interaction as before.

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


