***************
Atom selections
***************

Msys implements an atom selection language similar to that of VMD.
From Python, the language may be used to select atoms or clone subsets
of a ``System``, using the ``select``, ``selectIds``, or ``clone``
methods of ``System``.  In C++, the ``Atomselect`` function provides
all atom selection functionality.

Grammar
--------

The Msys atom selection grammar supports several primitive types which 
can be combined in various ways using logical operators.  Primitive
types include the following:

* Keyword selections: an attribute followed by one or more values, ranges,
  or regular expressions::

    name CA             # a single value
    resid 10 20 30      # multiple values
    resid 10 to 30      # a range of values
    name "C.*"          # regular expression, recognized by double quotes
    resid 5 8 to 10     # combination of single value and a range

* Singleword selections: A boolean attribute.  This includes both
  built-in singlewords such as ``protein``, as well as those defined
  as macros, such as ``acidic`` ::

    protein             # selects atoms identified as protein
    acidic              # expands to 'resname ASP GLU

* String functions: these are similar in form to keyword selections, but
  they use their arguments in special ways::

    smarts 'c[$(c[Ox1])]c'  # a smarts query.
    paramtype nonbonded HP  # query on the nonbonded type of a particle

* Comparisons: An inequality formed by two expressions, at least one of which
  should be a function of atom attributes:: 

    x > 0               # atoms in the positive x halfspace

From these primitive selections, more complex selections may be formed
using the following constructs.

* Boolean operators: AND, OR, NOT::

    protein and not hydrogen    # heavy protein atoms
    not oxygen and water        # same as "(not oxygen) and water"
    not oxygen or water         # same as "(not oxygen) or water"

* SAME ``attribute`` AS (``selection``)::

    same residue as name CA     # residues containing a CA atom

* ``Within`` selections of various flavors::

    within 3 of protein         # atoms within 3 of protein, including protein
    exwithin 3 of protein       # atoms within 3 of protein, but not protein
    pbwithin 3 of protein       # uses minimum image distance
    withinbonds 2 of resid 10   # atoms within two bonds of resid 10
     
* ``Nearest selections``::

    nearest 10 to residue 2     # nearest 10 atoms to any atom in residue 2
    pbnearest 10 to residue 2   # same, but with minimum image distance

In comparisons, expressions can be formed in the following ways.

* Numeric literals, and keywords with numeric types::

    x < 3                       # atoms whose x coordinate is less than 3

* Composition of expressions using various mathematical operators, 
  parentheses, and functions::

    x + y * z < 3               # the usual precedence rules apply
    sqr(x)/36 + sqr(z)/125 < 1  # an ellipsoidal cylinder


Differences with VMD
---------------------

Although the atom selection language in msys is similar to VMD's, there
are some important differences to bear in mind if you switch between
them:

* Element matching: In Msys, the atom selections "carbon", "hydrogen",
  "oxygen", etc. are based on the atomic number of the atoms.  In VMD,
  these atom selections are computed using a regular expression based
  on the atom name::

    vmd > atomselect macro oxygen name "O.*"

    vmd > atomselect macro hydrogen name "[0-9]?H.*"

    vmd > atomselect macro nitrogen name "N.*"

    vmd > atomselect macro carbon name "C.*" and not ion

* Implicit 'and': in VMD, selections can sometimes be concatenated with
  an implicit 'and'; e.g. "water within 3 of protein" will be parsed by
  VMD as "water and within 3 of protein".  In Msys, omitting the 'and' will
  result in a parse error.

* Field size: DMS and MAE files can hold chain, segment, and residue names
  of arbitrary length.  In Msys, these values are used as-is.  In VMD,
  the values are truncated; in particular, chain will be truncated to
  a single character in VMD, but not by Msys.

* Data representation: Msys has no concept of secondary structure, so the
  "sheet", helix", etc. atom selection keywords are not implemented in 
  msys.
  
* Floating-point roundoff: There may occasionally be differences in the
  results of distance based atom selections simply due the fact that Msys
  stores positions as doubles, while VMD stores them as floats.  

Built-in selections
*******************

The following selection keywords are available:

  ================  =========== ===========================================
  keyword           type        definition
  ================  =========== ===========================================
  atomicnumber      integer     `Atom`.atomic_number
  element           string      Abbreviation for element `Atom`.atomic_number
  chain             string      `Chain`.name
  segid             string      `Chain`.segid
  charge            float       `Atom`.charge
  fragment          integer     Connected residues will all have the same 
                                fragment id, except that the connection 
                                check will not follow disulfide bridges, 
                                identified as atoms whose name is "SG".
  index             integer     `Atom`.id
  mass              float       `Atom`.mass
  name              string      `Atom`.name
  numbonds          integer     `Atom`.nbonds - includes bonds to pseudoatoms
  degree            integer     number of bonds to real atoms; 0 for pseudoatoms
  resid             integer     `Residue`.resid
  residue           integer     `Residue`.id
  resname           string      `Residue`.name
  fragid            integer     `Atom`.fragid.  Connnected atoms will all 
                                have the same fragid.
  x, y, z           float       `Atom`.x, `Atom`.y, `Atom`.z, the position.
  vx, vy, vz        float       `Atom`.vx, `Atom`.vy, `Atom`.vz, the velocity.
  ================  =========== ===========================================


The following selection singlewords are available.  

  ===============   ==========================================================
  singleword        definition
  ===============   ==========================================================
  all               Every atom.
  none              No atoms.
  water             Atoms belonging to a residue containing the atomic number 
                    and bond structure of water, as well as those residues 
                    whose residue name is one of the following: "H2O", "HH0", 
                    "OHH", "HOH", "OH2", "SOL", "WAT", "TIP", "TIP2", "TIP3", 
                    "TIP4", "SPC".
  hydrogen           atomic number 1
  backbone          This singleword includes both protein backbone as well as 
                    nucleic backbone.  Protein backbone is identified by 
                    searching for atoms named "CA", "C", "O", and "N" in the 
                    same residue, and for atoms named "OT1", "OT2", "OXT", 
                    "O1", or "O2" that are bonded to one of the members of 
                    the first list.  If at least four such atoms are found, 
                    those atoms are identified as backbone.  Similarly, 
                    nucleic acid backbone atom names are P", "O1P", "O2P", 
                    "OP1", "OP2", "C3*", "C3'", "O3*", "O3'", "C4*", "C4'", 
                    "C5*", "C5'", "O5*", or "O5'"; or atoms named "H5T" or
                    "H3T" bonded to a member of the first set.  At least 
                    four such atoms must be found in the same residue in 
                    order to be identified as backbone.
  protein           residues containing protein backbone atoms.
  nucleic           residues containing nucleic backbone atoms.
  ===============   ==========================================================


The following are implemented as macros.

  ===========   ==========
  macro         definition
  ===========   ==========
  at            resname ADE A THY T
  acidic        resname ASP GLU
  cyclic        resname HIS PHE PRO TRP TYR
  acyclic       protein and not cyclic
  aliphatic     resname ALA GLY ILE LEU VAL
  alpha         protein and name CA
  amino         protein
  aromatic      resname HIS PHE TRP TYR
  basic         resname ARG HIS LYS HSP
  bonded        degree > 0
  buried        resname ALA LEU VAL ILE PHE CYS MET TRP
  cg            resname CYT C GUA G
  charged       basic or acidic
  hetero        not (protein or nucleic)
  hydrophobic   resname ALA LEU VAL ILE PRO PHE MET TRP
  small         resname ALA GLY SER
  medium        resname VAL THR ASP ASN PRO CYS ASX PCA HYP
  large         protein and not (small or medium)
  neutral       resname VAL PHE GLN TYR HIS CYS MET TRP ASX GLX PCA HYP
  polar         protein and not hydrophobic
  purine        resname ADE A GUA G
  pyrimidine    resname CYT C THY T URA U
  surface       protein and not buried
  lipid         resname DLPE DMPC DPPC GPC LPPC PALM PC PGCL POPC POPE
  lipids        lipid
  legacy_ion           resname AL BA CA Ca CAL CD CES CLA CL 'Cl-' Cl CO CS CU Cu CU1 CUA HG IN IOD K 'K+' MG MN3 MO3 MO4 MO5 MO6 NA Na NAW OC7 PB POT PT RB SOD TB TL WO4 YB ZN ZN1 ZN2
  ion           degree 0 and not atomicnumber 0 1 2 5 6 7 8 10 18 36 54 86
  ions          ion
  sugar         resname AGLC
  solvent       not (protein or sugar or nucleic or lipid)
  carbon        atomicnumber 6
  nitrogen      atomicnumber 7
  oxygen        atomicnumber 8
  sulfur        atomicnumber 16
  noh           not hydrogen
  heme          resname HEM HEME
  ===========   ==========


Smarts pattern selections
-------------------------
A SMARTS pattern is like a regular expression for molecular structures;
it's a concise way of specifying what sort of atom types and topology
you are looking for.  SMARTS patterns can be embedded in an atom selection
by providing the keyword 'smarts' followed by one or more SMARTS patterns,
which you will need to surround in single quotes if it contains any special
characters like parentheses::

    # select benzene rings
    mol.select("smarts 'c1ccccc1'")

See the description of the ``Smarts`` class for more information.


Parameter type selections
-------------------------

If a ParamTable contains a column named 'type', you can query for
atoms which participate in an interaction involving that type using the
'paramtype' keyword.  For example::

    # select atoms whose nonbonded type is 'H1'
    mol.select("paramtype nonbonded H1")


Comparison selections
---------------------

Comparisons are formed from two expressions and a binary comparison
operator.  The available comparison operators are the usual inequality
and equality operators: ``<``, ``>``, ``<=``, ``>=``, ``==``, and ``!=``.
Expressions can be built up from numeric literals and from keywords of 
float type, in the following ways:

* Binary mathematical operators: ``+``, ``-``, ``*``, and ``/``; e.g.,
  "x * y - z < 3".

* The C-style modulus function ``%``; e.g., "residue % 10 == 0" for every
  10th residue.

* Unary ``-``.

* The functions ``sqr``, ``sqrt``, and ``abs``; e.g., "sqrt(sqr(x)+sqr(y))<5".


User-defined keywords
---------------------

In addition to the aforementioned built-in keywords, any atom property may
also be used as an atom selection keyword.  For example::

  # add atom property 'foo' to a system.  The default value is empty string
  mol.addAtomProp('foo', str)

  # set the foo property to 'jrg' for all alpha carbons
  for a in mol.select('name CA'): a['foo'] = 'jrg'

  # check that selecting for foo equal to jrg is equivalent to 'name CA'
  assert mol.select('foo jrg') == mol.select('name CA')


User-defined atom selection macros
----------------------------------

This feature was removed in msys 1.7.7.

