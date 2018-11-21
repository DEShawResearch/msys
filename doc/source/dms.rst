
The DMS file format
===================

The preferred file format for Msys structures is the DMS file.


DMS versioning
--------------

Beginning with msys 1.7.0, a **dms_version** table is included in DMS
files written by msys.  The version table schema consists of a major
and minor version number, and will correspond to the major and minor
version of msys.  Going forward, msys will refuse to load DMS files
whose version is is higher than the msys version; thus, if and when
msys reaches version 1.8, files written by that version of msys will not
(necessarily) be readable by msys 1.7.  There is always the possibility
that forward compatibility could be ported to a later msys 1.7 version.
Backward compatibility with older dms versions will always be maintained.

The DMS versioning scheme serves to prevent problems arising from new
data structures being added to the DMS file in newer versions of msys 
which are not properly recognized by older versions.  For example,
the **nonbonded_combined_param** table was added in msys 1.4.0, but
because there was no dms version string at that time, older versions of
msys would have treated that table as an auxiliary table instead of
as a set of overrides to the nonbonded table.

In order to maintain this dms version semantics, it behooves the developers
of msys to think carefully about any new artifact appearing in the dms
file, or how data in the dms file is represented.  Any new dms table that
changes the forcefield or the result of atom selection must bring with
it an increment to the msys minor version.  

Schema changes in 1.7
^^^^^^^^^^^^^^^^^^^^^

An **msys_ct** table was added, which stores arbitrary key-value pairs
for components of the full system.  The **particle** table has a new
column called *ct* which assigns the particle to a row of the **msys_ct**
table.  

Msys 1.6.x should be able to read msys 1.7.x files with no problems, since
the ct information will just be ignored with no effect on atom selections
or forcefield.

DMS Schema
----------

All data in a DMS file lives in a flat list of two-dimensional tables.
Each table has a unique name.  Columns in the tables have a name, a
datatype, and several other attributes, most importantly, whether or
not the column is the primary key for the table.  Rows in the tables
hold a value for each of the columns.  Table names, column names, and
datatypes are case-preserving, but case-insensitive: thus 'pArTiCLE'
is the same table as 'particle', and 'NAME' is the same column as 'name'.

Of the five datatypes available in SQLite, DMS uses three: INTEGER, a
signed 64-bit int; FLOAT, a 64-bit IEEE floating point number; and TEXT,
a UTF8 string.  In addition, any value other than a primary key can be
NULL, indicating that no value is stored for that row and column.  A NULL
value is allowed in the DMS file but might be regarded as an invalid
value by a particular application; for example, Desmond make no use of
the atomic number column in the particle table, but Viparr requires it.

Because DMS is used to store dimensionful quantities, we must declare a
system of units.  The units in DMS reflect a compromise between an ideal
of full consistency and the reality of practical usage; in particular,
the mass unit is amu, rather than an algebraic combination of the energy,
length, and time units.

  =========     ========== 
  Dimension        Unit
  =========     ========== 
  TIME          picosecond 
  CHARGE        electron charge
  LENGTH        Angstrom
  ENERGY        thermochemical kcal/mol
  MASS          atomic mass unit (amu)
  =========     ========== 

In addition to tables, DMS files may contain stored queries known as views.
A view combines data from one or more tables, and may apply a predicate
as well a sorting criterion.  How this is actually done in SQL will be
obvious to database afficiandos; for this specification it suffices to
note that a view looks just like a table when reading a DMS file, so
the views will be described in terms of the data in their columns,
just as for tables.  Importantly, views cannot be written to directly;
one instead updates the tables to which they refer.

The DMS file contains the identity of all particles in the structure
as well as their positions and velocities in a global coordinate system.
The particle list includes both physical atoms as well as pseudoparticles
such as virtual sites and drude particles.  The most important table
has the name **particle**; all other tables containing additional particle
properties or particle-particle interactions refer back to records in
the **particle** table.  References to particles should follow a naming
convention of *p0*, *p1*, *p2*, ... for each particle referenced.

Particles
^^^^^^^^^

The **particle** table associates a unique *id* to all particles
in the structure.  The ids of the particles must all be contiguous,
starting at zero.  The ordering of the particles in a DMS file for the
purpose of, e.g., writing coordinate data, is given by the order of
their ids.  The minimal schema for the **particle** table is:

  ======    =======     ===========
  Column    Type        Description
  ======    =======     ===========
  anum      INTEGER     atomic number
  id        INTEGER     unique particle identifier   
  x         FLOAT       x-coordinate in LENGTH       
  y         FLOAT       y-coordinate in LENGTH       
  z         FLOAT       z-coordinate in LENGTH       
  ======    =======     ===========


Bonds
^^^^^

  ======    =======     ===========
  Column    Type        Description
  ======    =======     ===========
  p0        INTEGER     1st particle id 
  p1        INTEGER     2nd particle id 
  order     FLOAT       bond order      
  ======    =======     ===========

The **bond** table specifies the chemical topology of the system.  Here,
the topology is understood to be independent of the forcefield that describes
the interactions between particles.  Whether a water molecule is described
by a set of stretch and angle terms, or by a single constraint term, one would
still expect to find entries in the **bond** table corresponding to the
two oxygen-hydrogen bonds.  Bonds may also be present between a pseudoatom
and its parent particle or particles; these bonds aid in visualization.

The *p0* and *p1* values correspond to an id in the **particle** table.
Each *p0*, *p1* pair should be unique, non-NULL, and satisfy *p0 < p1*.

The global cell
^^^^^^^^^^^^^^^

  ======    =======     ===========
  Column    Type        Description
  ======    =======     ===========
  id        INTEGER     vector index (0, 1, or 2)    
  x         FLOAT       *x* component in LENGTH      
  y         FLOAT       *y* component in LENGTH      
  z         FLOAT       *z* component in LENGTH      
  ======    =======     ===========

The global_cell table specifies the dimensions of the periodic cell
in which particles interact.  There shall be three records, with *id*
0, 1, or 2; the primary key is provided since the order of the records
matters, and one would otherwise have difficulty referring to or updating
a particular record in the table.

Additional particle properties
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Additional per-particle properties not already specified in the
**particle** table should be added to the particle table as columns.
Here are the schema for the additional properties expected and/or
recognized by Desmond and by Viparr.

  ===============   =======     ===========
  Column            Type        Description
  ===============   =======     ===========
  mass              FLOAT       Desmond: particle mass in MASS        
  charge            FLOAT       Desmond: particle charge in CHARGE    
  vx                FLOAT       Desmond: x-velocity in LENGTH/TIME    
  vy                FLOAT       Desmond: y-velocity in LENGTH/TIME    
  vz                FLOAT       Desmond: z-velocity in LENGTH/TIME    
  nbtype            INTEGER     Desmond: nonbonded type 
  grp_temperature   INTEGER     Desmond: temperature group        
  grp_energy        INTEGER     Desmond: energy group             
  grp_ligand        INTEGER     Desmond: ligand group             
  grp_bias          INTEGER     Desmond: force biasing group      
  resid             INTEGER     Viparr: residue number               
  resname           TEXT        Viparr: residue name                 
  chain             TEXT        Viparr: chain identifier             
  name              TEXT        Viparr: atom name                    
  formal_charge     FLOAT       Viparr: format particle charge 
  occupancy         FLOAT       pdb occupancy value          
  bfactor           FLOAT       pdb temperature factor       
  ===============   =======     ===========

Forcefields
-----------

A description of a forcefield comprises the functional form of the
interactions between particles in a chemical system, the particles that
interact with a given functional form, and the parameters that govern a
particular interaction.  At a higher level, interactions can be described
as being \emph{local} or \emph{nonlocal}.  Local particle interactions in DMS
are those described by a fixed set of n-body terms.  These include bonded
terms, virtual sites, constraints, and polar terms.  Nonlocal interactions
in principle involve all particles in the system, though in practice
the potential is typically range-limited.  These include van der Waals
(vdw) interactions as well as electrostatics.  

Local particle interactions
^^^^^^^^^^^^^^^^^^^^^^^^^^^

In order to evaluate all the different forces between particles, a
program needs to be able to find them within a DMS file that may well
contain any number of other auxiliary tables.  The DMS format solves
this problem by providing a set of 'metatables' containing the names
of force terms required by the forcefield as well as the names of the
tables in which the force term data is found.  The force terms are placed
into one of four categories: bonded terms, constraints, virtual sites,
polar terms, described below.

  ===================   ===========
  Metatable name        Description
  ===================   ===========
  **bond_term**         Interactions representing bonds between atoms, including stretch, angle, and dihedral terms, as well as 1-4 pairs and position restraints.
  **constraint_term**   Constraints on bonds and/or angles involving a reduction in the number of degrees of freedom of the system.
  **virtual_term**      Similar to a constraint; a set of parameters describing how a pseudoparticle is to be positioned relative to a set of parent atoms. 
  **polar_term**        Similar to a virtual site; a set of parameters describing how a pseudoparticle moves relative to its parent atoms. 
  **nonbonded_table**   Additional or alternative nonbonded interactions.  Present only if such alternative tables are present.
  ===================   ===========

Each table name corresponding to the values in the local term metatables
is the designated string for a particular functional form.
The required columns for these tables is given in the next section.  Note
that creators of DMS files are free to implement the schema as an SQL
view, rather than as a pure table; a reader of a DMS file should not assume
anything about how the columns in the table name have been assembled.

Nonbonded interactions
^^^^^^^^^^^^^^^^^^^^^^

The functional form for nonbonded interactions, as well as the
tables containing the interaction parameters and type assignments,
are given by the fields in the **nonbonded_info** table, shown below:

  ===============   =======     ===========
  Column            Type        Description
  ===============   =======     ===========
  name              TEXT        nonbonded functional form 
  rule              TEXT        combining rule for nonbonded parameters 
  ===============   =======     ===========

There should exactly one record in the **nonbonded_info** table.
Like the local interaction tables,
the *name* field indicates the functional form of the nonbonded
interaction type.  If the particles have no nonbonded interactions,
*name* should have the special value `none`.

The parameters for nonbonded interactions will be stored in a table
called **nonbonded_param**, whose schema depends on the value of
*name* in **nonbonded_info**.  All such schemas must have a
primary key column called *id*; there are no other restrictions.

The *nbtype* column in the **particle** table gives the nonbonded
type assignment.  The value of the type assignment must correspond to
one of the primary keys in the **nonbonded_param** table.

Typically, the parameters governing the nonbonded interaction between
a pair of particles is a simple function of the parameters assigned to
the individual particles.  For example, in a Lennard-Jones functional
form with parameters *sigma* and *epsilon*, the combined parameters are
typically the arithmetic or geometric mean of *sigma* and *epsilon*.
The required approach is obtained by the application from the value of
*rule* in **nonbonded_info**.

For the interaction parameters that cannot be so simply derived, a table
called **nonbonded_combined_param** may be provided, with a schema shown
in Table~\ref{tab:combinedparam}.  Like the **nonbonded_param** table,
the schema of **nonbonded_combined_param** will depend on the functional
form of the nonbonded interactions, but there are two required columns,
which indicate which entry in **nonbonded_param** are being overridden.
Only *param1* and *param2* are required; the remaining columns provide
the interaction-dependent coefficients.

  ===============   =======     ===========
  Column            Type        Description
  ===============   =======     ===========
  param1            INTEGER      1st entry in **nonbonded_param** table
  param2            INTEGER      2nd entry in **nonbonded_param** table
  coeff1            FLOAT        first combined coefficient 
  \                              other combined coefficients... 
  ===============   =======     ===========


Alchemical systems
------------------

Methods for calculating relative free energies or energies of solvation
using free energy perturbation (FEP) involve mutating one or more chemical
entities from a reference state, labeled 'A', to a new state, labeled
'B'.  DMS treats FEP calculations as just another set of interactions
with an extended functional form.  In order to permit multiple independent
mutations to be carried out in the same simulation, a 'moiety' label is
applied to each mutating particle and bonded term.

Alchemical particles
^^^^^^^^^^^^^^^^^^^^

Any particle whose charge or nonbonded parameters changes in going
from state A to state B, is considered to be an alchemical particle
and must have a moiety assigned to it.  The set of distinct moieties
should begin at 0 and increase consecutively.  The set of alchemical
particles, if any, 
should be provided in a table called **alchemical_particle** shown
below:

  ===============   =======     ===========
  Column            Type        Description
  ===============   =======     ===========
  p0                INTEGER     alchemical particle id 
  moiety            INTEGER     moiety assignment 
  nbtypeA           INTEGER     entry in nonbonded_param for A state 
  nbtypeB           INTEGER     entry in nonbonded_param for B state 
  chargeA           FLOAT       charge in the A state 
  chargeB           FLOAT       charge in the B state 
  ===============   =======     ===========

Bonded terms
^^^^^^^^^^^^

Alchemical bonded terms are to be treated by creating a table analogous
to the non-alchemical version, but replacing each interaction parameter
with an 'A' and a 'B' version.  As a naming convention, the string
`alchemical_` should be prepended to the name of the table.  An example
is given below for **alchemical_stretch_harm** records, corresponding
to alchemical harmonic stretch terms with a functional form given by
interpolating between the parameters for states A and B.

  ===============   =======     ===========
  Column            Type        Description
  ===============   =======     ===========
  r0A               FLOAT       equilibrium separation in A state 
  fcA               FLOAT       force constant in A state 
  r0B               FLOAT       equilibrium separation in B state 
  fcB               FLOAT       force constant in B state 
  ---------------   -------     -----------
  p0                INTEGER     1st particle 
  p1                INTEGER     2nd particle 
  moiety            INTEGER     chemical group 
  ===============   =======     ===========

Constraint terms
^^^^^^^^^^^^^^^^

No support is offered for alchemical constraint terms at this time.
If particles A, b, and c are covered by an AH2 constraint in the A
state, and particles A, d, and e are covered by an AH2 constraint in
the B state, then the set of constraint terms in the alchemical DMS file
should include an AH4 constraint between A and b, c, d and e.

Virtual sites
^^^^^^^^^^^^^

No support is offered for alchemical virtual sites at this time.

Polar terms
^^^^^^^^^^^

No support is offered for alchemical polar terms at this time.


