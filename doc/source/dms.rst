*********
DMS Files
*********

Msys defines and implements an sqlite-based schema called DMS for
chemical systems.   This section provides an overview of the DMS
format which will be useful for users who wish to inspect their DMS
files manually using native sqlite tools, or who need to understand
the functional form of the forcefield tables found in DMS files in
order to, say, convert between file formats or use msys systems in
their own programs.

Overview
========

All data in a DMS file lives in a flat list of two-dimensional tables.
Each table has a unique name.  Columns in the tables have a name, a
datatype, and several other attributes, most importantly, whether or
not the column is the primary key for the table.  Rows in the tables
hold a value for each of the columns.  Table names, column names, and
datatypes are case-preserving, but case-insensitive: thus 'pArTiCLE'
is the same table as 'particle', and 'NAME' is the same column as 'name'.

In addition to tables, DMS files may contain stored queries known as views.
A view combines data from one or more tables, and may apply a predicate
as well a sorting criterion.  How this is actually done in SQL will be
obvious to database afficiandos; for this specification it suffices to
note that a view looks just like a table when reading a DMS file, so
the views will be described in terms of the data in their columns,
just as for tables.  Importantly, views cannot be written to directly;
one instead updates the tables to which they refer.

Units
-----

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

Versioning
----------

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


Chemical Structure
==================

The DMS file contains the identity of all particles in the structure
as well as their positions and velocities in a global coordinate system.
The particle list includes both physical atoms as well as pseudoparticles
such as virtual sites and drude particles.  The most important table
has the name **particle**; all other tables containing additional particle
properties or particle-particle interactions refer back to records in
the **particle** table.  References to particles should follow a naming
convention of *p0*, *p1*, *p2*, ... for each particle referenced.

Particles
---------

The **particle** table associates a unique *id* to all particles
in the structure.  The ids of the particles must all be contiguous,
starting at zero.  The ordering of the particles in a DMS file for the
purpose of, e.g., writing coordinate data, is given by the order of
their ids.  The minimal schema for the **particle** table is:

  ============= =======     ===========
  Column        Type        Description
  ============= =======     ===========
  anum          INTEGER     atomic number
  id            INTEGER     unique particle identifier   
  msys_ct       INTEGER     ct identifier
  x             FLOAT       x-coordinate in LENGTH       
  y             FLOAT       y-coordinate in LENGTH       
  z             FLOAT       z-coordinate in LENGTH       
  mass          FLOAT       particle mass in MASS        
  charge        FLOAT       particle charge in CHARGE    
  vx            FLOAT       x-velocity in LENGTH/TIME    
  vy            FLOAT       y-velocity in LENGTH/TIME    
  vz            FLOAT       z-velocity in LENGTH/TIME    
  nbtype        INTEGER     nonbonded type 
  resid         INTEGER     residue number               
  resname       TEXT        residue name                 
  chain         TEXT        chain identifier             
  name          TEXT        atom name                    
  formal_charge FLOAT       format particle charge 
  ============= =======     ===========


Msys organizes chemical system into a hierarchical structure.  The
hierarchy has the following names: ct, chain, residue, atom.  Rows
in the particle table of a dms file are mapped into these four
structural levels according to one or more columns in the particle
table.  The relevant columns for each structural level are:

    =========   =======
    structure   columns
    =========   =======
    ct          msys_ct
    chain       chain,segid
    residue     resname,resid,insertion
    atom        id
    =========   =======

Of these columns, only the id column is required.  The id will be
contiguous and start at 0.  The id determines the order of the particles
in the structure, important when dealing with simulation trajectories.
The other columns are treated as 0/empty string if not present.

Particles are mapped to ct object according to their msys_ct value.
Within a ct, there will be one chain object for each distinct
(chain,segid) tuple.  Within a chain object, there will be one residue
object for each distinct (resname,resid,insertion) tuple.  For example,
in the following hypothetical particle table with most columns elided:

    ==  =====   =====
    id  chain   resid   
    ==  =====   =====
    0   A       1
    1   A       1
    2   B       1
    3   C       2
    4   B       2
    ==  =====   =====

there would be one ct containing three chains with 1, 2, and 1 residues
in chains A, B and C, respectively.  Residues A/1, B/1, B/2, and C/2
would contain atoms 0-1, 2, 3 and 4.

Bonds
-----

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
---------------

  ======        =======     ===========
  Column        Type        Description
  ======        =======     ===========
  id            INTEGER     vector index (0, 1, or 2)    
  x             FLOAT       *x* component in LENGTH      
  y             FLOAT       *y* component in LENGTH      
  z             FLOAT       *z* component in LENGTH      
  ======        =======     ===========

The global_cell table specifies the dimensions of the periodic cell
in which particles interact.  There shall be three records, with *id*
0, 1, or 2; the primary key is provided since the order of the records
matters, and one would otherwise have difficulty referring to or updating
a particular record in the table.

Additional particle properties
------------------------------

Additional per-particle properties not already specified in the
**particle** table should be added to the particle table as columns.

  ===============   =======     ===========
  Column            Type        Description
  ===============   =======     ===========
  grp_temperature   INTEGER     temperature group        
  grp_energy        INTEGER     energy group             
  grp_ligand        INTEGER     ligand group             
  grp_bias          INTEGER     force biasing group      
  occupancy         FLOAT       pdb occupancy value          
  bfactor           FLOAT       pdb temperature factor       
  ===============   =======     ===========

Ct properties
-------------

The **msys_ct** table holds properties of each *ct* in the System.
The *msys_ct* field in the **particle** table maps each particle
to a ct.  The **msys_ct** table has only one required column,
*msys_name*, which holds the name of the ct.  Additional columns
are created in this table to hold ct properties.


Forcefields
===========

A description of a forcefield comprises the functional form of the
interactions between particles in a chemical system, the particles that
interact with a given functional form, and the parameters that govern a
particular interaction.  At a higher level, interactions can be described
as being `local` or `nonlocal`.  Local particle interactions in DMS
are those described by a fixed set of n-body terms.  These include bonded
terms, virtual sites, constraints, and polar terms.  Nonlocal interactions
in principle involve all particles in the system, though in practice
the potential is typically range-limited.  These include van der Waals
(vdw) interactions as well as electrostatics.  

Metatables
----------

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

Bond Terms
----------

Stretch terms
^^^^^^^^^^^^^

The vibrational motion between two atoms :math:`(i,j)` is represented by a
harmonic potential as:

.. math::

   V_s(r_{ij}) = f_c(r_{ij}-r_0)^2

where :math:`f_c` is the bond force constant in units of
:math:`\mathrm{Energy}/\mathrm{Length}^2` and :math:`r_0` is the equilibrium
bond distance.  Terms in ``stretch_harm``
are evaluated using this potential.

.. tabularcolumns:: |l|l|L|

.. _`tab:stretchharm`:

.. table:: Schema for the ``stretch_harm`` table

   +-----------+-------+----------------------------------------------+
   |name       |type   |description                                   |
   +===========+=======+==============================================+
   |r0         |FLOAT  |equilibrium separation (LENGTH)               |
   +-----------+-------+----------------------------------------------+
   |fc         |FLOAT  |force constant (ENERGY / LENGTH\ :sup:`2`)    |
   +-----------+-------+----------------------------------------------+
   |p0         |INTEGER|1st particle                                  |
   +-----------+-------+----------------------------------------------+
   |p1         |INTEGER|2nd particle                                  |
   +-----------+-------+----------------------------------------------+
   |constrained|INTEGER|if nonzero, constrained; default 0            |
   +-----------+-------+----------------------------------------------+

Stretch terms that overlap with
constraints should have the constrained field set to 1.  Applications
that evaluate constraint terms need not evaluate ``stretch_harm``
records that are marked as constrained.

Angle terms
^^^^^^^^^^^

The angle vibration between three atoms :math:`(i,j,k)`  is evaluated as:

.. math::

   V_a(\theta_{ijk}) = f_c(\theta_{ijk}-\theta_0)^2

where :math:`f_c` is the angle force constant in
:math:`\mathrm{Energy}/\mathrm{Radians}^2` and :math:`\theta_0` is the equilibrium
angle in radians.  Beware, the explicit use of the :math:`\theta_{ijk}`
angle will introduce discontinuities in the potential at
:math:`\theta_{ijk} = \pm\pi`\ .  Terms in ``angle_harm``
are evaluated using this potential.

.. tabularcolumns:: |l|l|L|

.. `tab:angleharm`::

.. table:: Schema for the ``angle_harm`` table

   +-----------+-------+----------------------------------------------+
   |name       |type   |description                                   |
   +===========+=======+==============================================+
   |theta0     |FLOAT  |equilibrium angle (DEGREES)                   |
   +-----------+-------+----------------------------------------------+
   |fc         |FLOAT  |force constant (ENERGY / RADIAN\ :sup:`2`)    |
   +-----------+-------+----------------------------------------------+
   |p0         |INTEGER|1st particle                                  |
   +-----------+-------+----------------------------------------------+
   |p1         |INTEGER|2nd particle                                  |
   +-----------+-------+----------------------------------------------+
   |p2         |INTEGER|3rd particle                                  |
   +-----------+-------+----------------------------------------------+
   |constrained|INTEGER|constrained if nonzero; default 0             |
   +-----------+-------+----------------------------------------------+

The :math:`p0` particle forms the
vertex.  Angle terms that overlap with constraints should have the
constrained field set to 1.  Applications that evaluate constraint
terms need not evaluate ``angle_harm`` records that are marked as
constrained.

Proper dihedral terms
^^^^^^^^^^^^^^^^^^^^^

Two functional forms for calculating proper and improper torsion
potential terms are specified.  The first is:

.. math::

   V_t(\phi_{ijkl}) = f_{c0}
   + \sum_{n=1}^6 f_{cn} \cos(n\phi_{ijkl}-\phi_0)

where :math:`f_{c0} \ldots f_{c6}` are dihedral angle force constants in units
of Energy and :math:`\phi_0` is the equilibrium dihedral angle
in radians.
The :math:`\phi` angle is formed by the planes :math:`p0`\ --\ :math:`p1`\ --\ :math:`p2` and
:math:`p1`\ --\ :math:`p2`\ --\ :math:`p3`\ .
Terms in ``dihedral_trig`` are handled by this potential function.

.. tabularcolumns:: |l|l|L|

.. _`tab:dihedral_trig`:

.. table:: Schema for the ``dihedral_trig`` table.

   +----+-------+-------------------------------+
   |name|type   |description                    |
   +====+=======+===============================+
   |phi0|FLOAT  |phase (DEGREES)                |
   +----+-------+-------------------------------+
   |fc0 |FLOAT  |order-0 force constant (ENERGY)|
   +----+-------+-------------------------------+
   |fc1 |FLOAT  |order-1 force constant (ENERGY)|
   +----+-------+-------------------------------+
   |fc2 |FLOAT  |order-2 force constant (ENERGY)|
   +----+-------+-------------------------------+
   |fc3 |FLOAT  |order-3 force constant (ENERGY)|
   +----+-------+-------------------------------+
   |fc4 |FLOAT  |order-4 force constant (ENERGY)|
   +----+-------+-------------------------------+
   |fc5 |FLOAT  |order-5 force constant (ENERGY)|
   +----+-------+-------------------------------+
   |fc6 |FLOAT  |order-6 force constant (ENERGY)|
   +----+-------+-------------------------------+
   |p0  |INTEGER|1st particle                   |
   +----+-------+-------------------------------+
   |p1  |INTEGER|2nd particle                   |
   +----+-------+-------------------------------+
   |p2  |INTEGER|3rd particle                   |
   +----+-------+-------------------------------+
   |p3  |INTEGER|4th particle                   |
   +----+-------+-------------------------------+

Improper dihedral terms
^^^^^^^^^^^^^^^^^^^^^^^

The second dihedral functional form is:

.. math::
   :label:  eqn:improper_harm

   V_t(\phi_{ijkl}) = f_c (\phi_{ijkl}-\phi_0)^2

where :math:`f_c` is the dihedral angle force constant in units of
Energy/radians\ :math:`^2` and :math:`\phi_0` is the equilibrium dihedral angle
in radians.  The :math:`\phi` angle is formed by the planes
:math:`p0`\ --\ :math:`p1`\ --\ :math:`p2` and :math:`p1`\ --\ :math:`p2`\ --\ :math:`p3`\ .  Terms in
``improper_harm`` are handled by this potential function.

The harmonic dihedral term given in Equation :eq:`eqn:improper_harm` can lead to
accuracy issues if :math:`f_c` is too small, or if initial conditions are poorly
chosen due to a discontinuity in the definition of the first derivative
with respect to :math:`i` in :math:`\phi_{ijkl}` near :math:`\phi_0 \pm \pi`\ .

.. tabularcolumns:: |l|l|L|

.. _`tab:improper_harm`:

.. table:: Schema for the ``improper_harm`` table.

   +----+-------+------------------------------------------+
   |name|type   |description                               |
   +====+=======+==========================================+
   |phi0|FLOAT  |equilibrium separation (DEGREES)          |
   +----+-------+------------------------------------------+
   |fc  |FLOAT  |force constant (ENERGY / DEGREE\ :sup:`2`)|
   +----+-------+------------------------------------------+
   |p0  |INTEGER|1st particle                              |
   +----+-------+------------------------------------------+
   |p1  |INTEGER|2nd particle                              |
   +----+-------+------------------------------------------+
   |p2  |INTEGER|3rd particle                              |
   +----+-------+------------------------------------------+
   |p3  |INTEGER|4th particle                              |
   +----+-------+------------------------------------------+

CMAP torsion terms
^^^^^^^^^^^^^^^^^^


CMAP is a torsion-torsion cross-term that uses a tabulated energy
correction.  It is found in more recent versions of the CHARMM
forcefield.  The potential function is given by:

.. math::

   V_c(\phi,\psi) =
   \sum_{n=1}^4\sum_{m=1}^4
   C_{nm}\left(\frac{\psi-\psi_L}{\Delta_\psi}\right)^{n-1}\left(\frac{\phi-\phi_L}{\Delta_\phi}\right)^{m-1}

where :math:`C_{nm}` are bi-cubic interpolation coefficients derived from
the supplied energy table, :math:`\phi` is the dihedral angle formed by
particles :math:`p0 \ldots p3`\ , and :math:`\psi` is the dihedral angle formed by particles
:math:`p4 \ldots p7`\ .  The grid spacings are also derived from
the supplied energy table.  Terms in ``torsiontorsion_cmap``
are handled by this potential function.

The ``cmap`` tables for each term can be found in ``cmapN``\ ,
where ``N`` is a unique integer identifier for a particular table
(multiple ``cmap`` terms in ``torsiontorsion_cmap`` can refer to a
single ``cmapN`` block).  The format of the ``cmap`` tables
consists of two torsion angles in degrees and an associated energy.
``cmap`` tables must begin with both torsion angles equal to -180.0 and
increase fastest in the second torsion angle.  The grid spacing must
be uniform within each torsion coordinate, but can be different from
the grid spacing in other torsion coordinates.
More information can
be found in [Bro-2004]_.

.. tabularcolumns:: |l|l|L|

.. _`tab:cmaptable`:

.. table:: Schema for each of the tables holding the 2D ``cmap`` grids

   +------+-----+------------------------+
   |name  |type |description             |
   +======+=====+========================+
   |phi   |FLOAT|phi coordinate (DEGREES)|
   +------+-----+------------------------+
   |psi   |FLOAT|psi coordinate (DEGREES)|
   +------+-----+------------------------+
   |energy|FLOAT|energy value (ENERGY)   |
   +------+-----+------------------------+

The CHARMM27 forcefield uses six ``cmap`` tables, which have names
``cmap1``\ , ``cmap2``\ , ..., ``cmap6`` in DMS.

.. tabularcolumns:: |l|l|L|

.. _`tab:torsiontorsion_cmap`:

.. table:: Schema for the ``torsiontorsion_cmap`` table

   +----+-------+------------------+
   |name|type   |description       |
   +====+=======+==================+
   |cmap|INTEGER|name of cmap table|
   +----+-------+------------------+
   |p0  |INTEGER|1st particle      |
   +----+-------+------------------+
   |p1  |INTEGER|2nd particle      |
   +----+-------+------------------+
   |p2  |INTEGER|3rd particle      |
   +----+-------+------------------+
   |p3  |INTEGER|4th particle      |
   +----+-------+------------------+
   |p4  |INTEGER|5th particle      |
   +----+-------+------------------+
   |p5  |INTEGER|6th particle      |
   +----+-------+------------------+
   |p6  |INTEGER|7th particle      |
   +----+-------+------------------+
   |p7  |INTEGER|8th particle      |
   +----+-------+------------------+

.. _`sec:position_restraints`:

Position restraint terms
^^^^^^^^^^^^^^^^^^^^^^^^

Particles can be restrained to a given global coordinate by means of
the restraining potential:

.. math::

   V_r(x,y,z) = \frac{\lambda}{2} (
   f_{cx}(x-x_0)^2
   + f_{cy}(y-y_0)^2
   + f_{cz}(z-z_0)^2
   )

where :math:`f_{cx}`\ , :math:`f_{cy}`\ , :math:`f_{cz}` are the force constants in
:math:`\mathrm{Energy}/\mathrm{Length}\sp{2}` and :math:`x_0`\ , :math:`y_0`\ , :math:`z_0` are the desired global
cell coordinates (units of Length).  :math:`\lambda` is a pure scaling factor, set to 1 by
default.  Terms in ``posre_harm`` are evaluated using this potential.

.. tabularcolumns:: |l|l|L|

.. _`tab:posre_harm`:

.. table:: Schema for the ``posre_harm`` table

   +----+-------+---------------------------------------------+
   |name|type   |description                                  |
   +====+=======+=============================================+
   |fcx |FLOAT  |X force constant in ENERGY/LENGTH\ :sup:`2`  |
   +----+-------+---------------------------------------------+
   |fcy |FLOAT  |Y force constant in ENERGY/LENGTH\ :sup:`2`  |
   +----+-------+---------------------------------------------+
   |fcz |FLOAT  |Z force constant in ENERGY/LENGTH\ :sup:`2`  |
   +----+-------+---------------------------------------------+
   |p0  |INTEGER|restrained particle                          |
   +----+-------+---------------------------------------------+
   |x0  |FLOAT  |x reference coordinate                       |
   +----+-------+---------------------------------------------+
   |y0  |FLOAT  |y reference coordinate                       |
   +----+-------+---------------------------------------------+
   |z0  |FLOAT  |z reference coordinate                       |
   +----+-------+---------------------------------------------+


Pair 12--6 terms
^^^^^^^^^^^^^^^^

Pair terms in ``pair_12_6_es`` allow for modifying the normally
calculated nonbonded interactions either by scaling the interaction
energy, or by specifying new coefficients to use for a particular
pair.  This partial or modified energy is calculated in addition to
the normally calculated interaction energy.

The functional form of the pair potential is:

.. math::

   V_p(r_{ij}) =
   \frac{a_{ij}}{r_{ij}^{12}}
   - \frac{b_{ij}}{r_{ij}^{ 6}}
   + \frac{q_{ij}}{r_{ij}}

The  :math:`a_{ij}`\ , :math:`b_{ij}`\ , and :math:`q_{ij}`  coefficients are specified
in the ``pair_12_6_es`` table.

.. tabularcolumns:: |l|l|L|

.. _`tab:pair_12_6_es`:

.. table:: Schema for the ``pair_12_6_es`` table

   +----+-------+-------------------------------------------------+
   |name|type   |description                                      |
   +====+=======+=================================================+
   |aij |FLOAT  |scaled LJ12 coeff in ENERGY LENGTH\ :sup:`12`    |
   +----+-------+-------------------------------------------------+
   |bij |FLOAT  |scaled LJ6 coeff in ENERGY LENGTH\ :sup:`6`      |
   +----+-------+-------------------------------------------------+
   |qij |FLOAT  |scaled product of charges in CHARGE\ :sup:`2`    |
   +----+-------+-------------------------------------------------+
   |p0  |INTEGER|1st particle                                     |
   +----+-------+-------------------------------------------------+
   |p1  |INTEGER|2nd particle                                     |
   +----+-------+-------------------------------------------------+

Flat-bottomed harmonic well
^^^^^^^^^^^^^^^^^^^^^^^^^^^

The functional form of the flat-bottomed harmonic angle term is
:math:`V=|d|^2` where

.. math::
    
    d &= \begin{cases}
    (\theta-\theta_0+\sigma) & \mbox{where } \theta-\theta_0 < -\sigma \\
    0                   & \mbox{where } -\sigma <= \theta-\theta_0 < \sigma \\
    (\theta-\theta_0-\sigma) & \mbox{where }  \sigma <= \theta-\theta_0 
    \end{cases}

and :math:`\theta_0` is in radians.

.. tabularcolumns:: |l|l|L|

.. _`tab:angle_fbhw`:

.. table:: Schema for the ``angle_fbhw`` table

   +------+-------+---------------------------------------------+
   |name  |type   |description                                  |
   +======+=======+=============================================+
   |fc    |FLOAT  |force constant in ENERGY/RADIANS\ :sup:`2`   |
   +------+-------+---------------------------------------------+
   |theta0|FLOAT  |equilibrium angle in DEGREES                 |
   +------+-------+---------------------------------------------+
   |sigma |FLOAT  |half-width of flat-bottomed region in DEGREES|
   +------+-------+---------------------------------------------+
   |p0    |INTEGER|first particle                               |
   +------+-------+---------------------------------------------+
   |p1    |INTEGER|second particle                              |
   +------+-------+---------------------------------------------+
   |p2    |INTEGER|third particle                               |
   +------+-------+---------------------------------------------+


The functional form of the FBHW improper term is :math:`V=f_c d^2` where

.. math::

    d &= \begin{cases}
    (\phi-\phi_0+\sigma) & \mbox{where } \phi-\phi_0 < -\sigma \\
    0                   & \mbox{where } -\sigma <= \phi-\phi_0 < \sigma \\
    (\phi-\phi_0-\sigma) & \mbox{where }  \sigma <= \phi-\phi_0 
    \end{cases}


The improper dihedral angle phi is the angle between the plane
ijk and jkl. Thus fc is in ENERGY and phi0 is in RADIANS. 


.. tabularcolumns:: |l|l|L|

.. _`tab:improper_fbhw`:

.. table:: Schema for the ``improper_fbhw`` table

   +-----+-------+----------------------------------------------+
   |name |type   |description                                   |
   +=====+=======+==============================================+
   |fc   |FLOAT  |force constant in ENERGY/RADIANS\ :sup:`2`    |
   +-----+-------+----------------------------------------------+
   |phi0 |FLOAT  |equilibrium improper dihedral angle in DEGREES|
   +-----+-------+----------------------------------------------+
   |sigma|FLOAT  |half-width of flat-bottomed region in DEGREES |
   +-----+-------+----------------------------------------------+
   |p0   |INTEGER|first particle                                |
   +-----+-------+----------------------------------------------+
   |p1   |INTEGER|second particle                               |
   +-----+-------+----------------------------------------------+
   |p2   |INTEGER|third particle                                |
   +-----+-------+----------------------------------------------+
   |p3   |INTEGER|fourth particle                               |
   +-----+-------+----------------------------------------------+


The functional form of the FBHW posre term is :math:`V=f_c/2 d^2` where

.. math::

    d = \begin{cases}
    |r-r0|-\sigma & \mbox{where } |r-r0| > \sigma \\
        0         & \mbox{where } |r-r0| <= \sigma
    \end{cases}


This is not as general as the fully harmonic position restraint term
in that you can't specify different force constants for the three
coordinate axes.

.. tabularcolumns:: |l|l|L|

.. _`tab:posre_fbhw`:

.. table:: Schema for the ``posre_fbhw`` table

   +-----+-------+-------------------------------------------+
   |name |type   |description                                |
   +=====+=======+===========================================+
   |fc   |FLOAT  |force constant in ENERGY/LENGTH\ :sup:`2`  |
   +-----+-------+-------------------------------------------+
   |x0   |FLOAT  |equilibrium :math:`x` coordinate in LENGTH |
   +-----+-------+-------------------------------------------+
   |y0   |FLOAT  |equilibrium :math:`y` coordinate in LENGTH |
   +-----+-------+-------------------------------------------+
   |z0   |FLOAT  |equilibrium :math:`z` coordinate in LENGTH |
   +-----+-------+-------------------------------------------+
   |sigma|FLOAT  |radius of flat-bottomed region in LENGTH   |
   +-----+-------+-------------------------------------------+
   |p0   |INTEGER|restrained particle                        |
   +-----+-------+-------------------------------------------+

Exclusions
----------

Exclusion terms in ``exclusion`` are used to prevent
calculation of certain non bonded interactions at short ranges.  The
excluded interactions are typically those that involve particles
separated by one or two bonds, as these interactions are assumed to be
adequately modeled by the stretch and angle terms described above.

.. _`tab:exclusion`:

.. tabularcolumns:: |l|l|L|

.. table:: Schema for the ``exclusion`` table

   +----+-------+------------+
   |name|type   |description |
   +====+=======+============+
   |p0  |INTEGER|1st particle|
   +----+-------+------------+
   |p1  |INTEGER|2nd particle|
   +----+-------+------------+

It is required that :math:`p0 < p1` for each term, and every :math:`p0`\ , :math:`p1`
pair should be unique.

Constraint Terms
----------------

Constraints fix the distances between pairs of particles according to
a topology of rigid rods:

.. math::

   || r_i - r_j || &= d_{ij}

   || r_k - r_l || &= d_{kl}

   \ldots

The topologies that can be constrained are:

+ **AHn**: n particles connected to a single particle, with
  :math:`1 \le n \le 8`.
+ **HOH**: three mutually connected particles.

The schemas in the DMS file for ``AHn`` and ``HOH`` constraints
are shown in :ref:`tab:constraint-ahn` and :ref:`tab:constraint-hoh`,
respectively.  Each record in the ``AHn`` table gives the length of
the bonds between a single parent atom and ``n`` child atoms.
Each record in the ``HOH`` table gives the angle between the two
O-H bonds and the respective bonds lengths.

.. tabularcolumns:: |l|l|L|

.. _`tab:constraint-ahn`:

.. table:: Schema for the ``constraint_ahN`` tables


   +------+-------+-----------------+
   |name  |type   |description      |
   +======+=======+=================+
   |r1    |FLOAT  |A-H1 distance    |
   +------+-------+-----------------+
   |r2    |FLOAT  |A-H2 distance    |
   +------+-------+-----------------+
   |...   |       |                 |
   +------+-------+-----------------+
   |rN    |FLOAT  |A-HN distance    |
   +------+-------+-----------------+
   |p0    |INTEGER|id of parent atom|
   +------+-------+-----------------+
   |p1    |INTEGER|id of H1         |
   +------+-------+-----------------+
   |p2    |INTEGER|id of H2         |
   +------+-------+-----------------+
   |...   |       |                 |
   +------+-------+-----------------+
   |pN    |INTEGER|id of HN         |
   +------+-------+-----------------+

.. tabularcolumns:: |l|l|L|

.. _`tab:constraint-hoh`:

.. table:: Schema for the ``constraint_hoh`` (rigid water) table

   +-----+-------+-------------------------+
   |name |type   |description              |
   +=====+=======+=========================+
   |theta|FLOAT  |H-O-H angle in DEGREES   |
   +-----+-------+-------------------------+
   |r1   |FLOAT  |O-H1 distance            |
   +-----+-------+-------------------------+
   |r2   |FLOAT  |O-H2 distance            |
   +-----+-------+-------------------------+
   |p0   |INTEGER|id of heavy atom (oxygen)|
   +-----+-------+-------------------------+
   |p1   |INTEGER|id of H1                 |
   +-----+-------+-------------------------+
   |p2   |INTEGER|id of H2                 |
   +-----+-------+-------------------------+

A constrained particle is no longer free; each such particle has
:math:`3 - m/2` degrees of freedom, where :math:`m` is the number of
independent constraints involved; for example, a pair of particles
having only one distance constraint between them has five degrees of
freedom.  Constraints thus affect the calculation of the instantaneous
temperature and pressure, which depend on the number of degrees of
freedom.

The ``AHnR`` constraints are versions of the ``AHn`` constraints with
additional distances sufficient to create a rigid body.  As with water
constraints, the alternative Reich algorithm is used by default.

.. tabularcolumns:: |l|l|L|

.. _`tab:constraint-ah1r`:

.. table:: Schema for the ``constraint_ah1R`` table

   +------+-------+-----------------+
   |name  |type   |description      |
   +======+=======+=================+
   |r1    |FLOAT  |A-H1 distance    |
   +------+-------+-----------------+
   |p0    |INTEGER|id of parent atom|
   +------+-------+-----------------+
   |p1    |INTEGER|id of H1         |
   +------+-------+-----------------+

.. tabularcolumns:: |l|l|L|

.. _`tab:constraint-ah2r`:


.. table:: Schema for the ``constraint_ah2R`` table

   +------+-------+-----------------+
   |name  |type   |description      |
   +======+=======+=================+
   |r1    |FLOAT  |A-H1 distance    |
   +------+-------+-----------------+
   |r2    |FLOAT  |A-H2 distance    |
   +------+-------+-----------------+
   |r3    |FLOAT  |H1-H2 distance   |
   +------+-------+-----------------+
   |p0    |INTEGER|id of parent atom|
   +------+-------+-----------------+
   |p1    |INTEGER|id of H1         |
   +------+-------+-----------------+
   |p2    |INTEGER|id of H2         |
   +------+-------+-----------------+

.. tabularcolumns:: |l|l|L|

.. _`tab:constraint-ah3r`:

.. table:: Schema for the ``constraint_ah3R`` table

   +------+-------+-----------------+
   |name  |type   |description      |
   +======+=======+=================+
   |r1    |FLOAT  |A-H1 distance    |
   +------+-------+-----------------+
   |r2    |FLOAT  |A-H2 distance    |
   +------+-------+-----------------+
   |r3    |FLOAT  |A-H3 distance    |
   +------+-------+-----------------+
   |r4    |FLOAT  |H1-H2 distance   |
   +------+-------+-----------------+
   |r5    |FLOAT  |H1-H3 distance   |
   +------+-------+-----------------+
   |r6    |FLOAT  |H2-H3 distance   |
   +------+-------+-----------------+
   |p0    |INTEGER|id of parent atom|
   +------+-------+-----------------+
   |p1    |INTEGER|id of H1         |
   +------+-------+-----------------+
   |p2    |INTEGER|id of H2         |
   +------+-------+-----------------+
   |p3    |INTEGER|id of H3         |
   +------+-------+-----------------+


Virtual sites
-------------

.. parsed-literal::

   force.virtual = {
     exclude = [ ... ] # optional names to remove
     include = [ ... ] # optional names to add
     # typically, no other per-term arguments required
   }

Virtual sites, a form of pseudoparticle, are additional off-atom
interaction sites that can be added to a molecular system.  These
sites can have charge or van der Waals parameters associated with
them; they are usually massless.  The TIP4P and TIP5P water models are
examples that contain one and two off-atom (virtual) sites,
respectively.  Because these sites are massless, it is necessary to
redistribute any forces acting on them to the particles used in their
construction.  (A consistent way to do this can be found in [Gun-1984]_.) The
virial in most cases must also be modified after redistributing the
virtual site force.

The types of virtual site placement routines are described below.

lc2 virtual site
^^^^^^^^^^^^^^^^

The lc2 virtual site is placed some fraction a along the vector
between two particles :math:`(i,j)`\ .

.. math::

   \vec r_v = (1-c_1)\vec r_i + c_1 \vec r_j

.. tabularcolumns:: |l|l|L|

.. _`tab:virtuallc2`:

.. table:: Schema for ``virtual_lc2`` records

   +----+-------+-----------------+
   |name|type   |description      |
   +====+=======+=================+
   |c1  |FLOAT  |coefficient 1    |
   +----+-------+-----------------+
   |p0  |INTEGER|pseudoparticle id|
   +----+-------+-----------------+
   |p1  |INTEGER|parent atom i    |
   +----+-------+-----------------+
   |p2  |INTEGER|parent atom j    |
   +----+-------+-----------------+

Pseudoparticle :math:`p0` is
placed at the fractional position :math:`c1` along the interpolated line between
:math:`p1` and :math:`p2`\ .

lc3 virtual site
^^^^^^^^^^^^^^^^

The lc3 virtual site is placed some fraction :math:`a` and :math:`b` along the
vectors between particles :math:`(i,j)` and :math:`(i,k)` respectively.  The
virtual particle lies in the plane formed by :math:`(i,j,k)`\ .

.. math::

   \vec r_v = (1-c_1-c_2)\vec r_i + c_1 \vec r_j + c_2 \vec r_k

.. tabularcolumns:: |l|l|L|

.. _`tab:virtual_lc3`:

.. table:: Schema for the ``virtual_lc3`` table

   +----+-------+-----------------+
   |name|type   |description      |
   +====+=======+=================+
   |c1  |FLOAT  |coefficient 1    |
   +----+-------+-----------------+
   |c2  |FLOAT  |coefficient 2    |
   +----+-------+-----------------+
   |p0  |INTEGER|pseudoparticle id|
   +----+-------+-----------------+
   |p1  |INTEGER|parent atom i    |
   +----+-------+-----------------+
   |p2  |INTEGER|parent atom j    |
   +----+-------+-----------------+
   |p3  |INTEGER|parent atom k    |
   +----+-------+-----------------+

fdat3 virtual site
^^^^^^^^^^^^^^^^^^

The fdat3 virtual site is placed at a fixed distance  :math:`d`  from
particle  :math:`i`\ , at a fixed angle :math:`\theta`
defined by particles  :math:`(v,i,j)`  and at a fixed torsion :math:`\phi`
defined by particles :math:`(v,i,j,k)`\ .

.. math::

   \vec r_v = \vec r_i +
   a \vec r_{1} +  b \vec r_{2} + c \vec r_{2}\times \vec r_{1}

where :math:`\vec r_1` and :math:`\vec r_2` are unit vectors defined by

.. math::

   \vec r_1 &\propto \vec r_j - \vec r_i

   \vec r_2 &\propto \vec r_k - \vec r_j - (\vec r_k - \vec r_j)\cdot \vec r_1 \vec r_1

The coefficients :math:`a`\ , :math:`b` and :math:`c` above are defined as
:math:`a = d\cos(\theta)`\ , :math:`b = d\sin(\theta)\cos(\phi)` and
:math:`c = d\sin(\theta)\sin(\phi)`\ .

.. tabularcolumns:: |l|l|L|

.. _`tab:virtual_fdat3`:

.. table:: Schema for the ``virtual_fdat3`` table

   +----+-------+--------------------------+
   |name|type   |description               |
   +====+=======+==========================+
   |c1  |FLOAT  |:math:`d` coefficient     |
   +----+-------+--------------------------+
   |c2  |FLOAT  |:math:`\theta` coefficient|
   +----+-------+--------------------------+
   |c3  |FLOAT  |:math:`\phi` coefficient  |
   +----+-------+--------------------------+
   |p0  |INTEGER|pseudoparticle id         |
   +----+-------+--------------------------+
   |p1  |INTEGER|parent atom i             |
   +----+-------+--------------------------+
   |p2  |INTEGER|parent atom j             |
   +----+-------+--------------------------+
   |p3  |INTEGER|parent atom k             |
   +----+-------+--------------------------+

out3 virtual site
^^^^^^^^^^^^^^^^^

The out3 virtual site can be placed out of the plane of three
particles :math:`(i,j,k)`\ .

.. math::

   \vec r_v = \vec r_i + c_1 (\vec r_j-\vec r_i) + c_2 (\vec r_k-\vec r_i)
   + c_3 (\vec r_j-\vec r_i) \times (\vec r_k-\vec r_i)

.. tabularcolumns:: |l|l|L|

.. _`tab:virtual_out3`:

.. table:: Schema for the ``virtual_out3`` table

   +----+-------+-----------------+
   |name|type   |description      |
   +====+=======+=================+
   |c1  |FLOAT  |coefficient 1    |
   +----+-------+-----------------+
   |c2  |FLOAT  |coefficient 2    |
   +----+-------+-----------------+
   |c3  |FLOAT  |coefficient 3    |
   +----+-------+-----------------+
   |p0  |INTEGER|pseudoparticle id|
   +----+-------+-----------------+
   |p1  |INTEGER|parent atom i    |
   +----+-------+-----------------+
   |p2  |INTEGER|parent atom j    |
   +----+-------+-----------------+
   |p3  |INTEGER|parent atom k    |
   +----+-------+-----------------+

Nonbonded interactions
----------------------

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

.. tabularcolumns:: |l|l|L|

.. _`tab:vdw126`:

.. table:: Schema for the **vdw_12_6** nonbonded type

   +-------+-----+--------------------+
   |name   |type |description         |
   +=======+=====+====================+
   |sigma  |FLOAT|VdW radius in LENGTH|
   +-------+-----+--------------------+
   |epsilon|FLOAT|VdW energy in ENERGY|
   +-------+-----+--------------------+

The functional form is
:math:`V = a_{ij} / |r|^{12} + b_{ij} / |r|^6`\ , where :math:`a_{ij}` and :math:`b_{ij}` are computed
by applying either the combining rule from ``nonbonded_info`` or the
value from ``nonbonded_combined_param`` to obtain
:math:`\sigma` and :math:`\epsilon`\ , then computing
:math:`a_{ij} = 4 \epsilon \sigma^{12}` and :math:`b_{ij} = -4 \epsilon \sigma^6`\ .



Alchemical systems
==================

Methods for calculating relative free energies or energies of solvation
using free energy perturbation (FEP) involve mutating one or more chemical
entities from a reference state, labeled 'A', to a new state, labeled
'B'.  DMS treats FEP calculations as just another set of interactions
with an extended functional form.  In order to permit multiple independent
mutations to be carried out in the same simulation, a 'moiety' label is
applied to each mutating particle and bonded term.

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

Alchemical bonded terms can be treated by creating a table analogous
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

References
==========

.. [Bro-2004] C. L. Brooks III, A. D. MacKerell Jr., M. Feig,
   "Extending the treatment of backbone energetics in protein force
   fields: limitations of gas-phase quantum mechanics in reproducing
   protein conformational distributions in molecular dynamics
   simulations", *J. Comput. Chem.*, 25:1400--1415, 2004.

.. [Gun-1984] W. F. van Gunsteren H. J. C Berendsen, "Molecular
   dynamics simulations: Techniques and approaches", In A. J. Barnes
   et al., editor, *Molecular Liquids: Dynamics and Interactions*,
   NATO ASI C 135, pages 475--500. Reidel Dordrecht, The Netherlands,
   1984.

