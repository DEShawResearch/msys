***************
Atom selections
***************

Msys implements essentially all of the atom selection language of VMD.
Differences between Msys and VMD's implementations (other than as yet
undiscovered bugs in Msys!) fall into the following categories:

* Element matching: In Msys, the atom selections "carbon", "hydrogen",
  "oxygen", etc. are based on the atomic number of the atoms.  In 
  VMD, you maybe shocked and surprised to learn, these atom selections
  are computed using a regular expression based on the atom name::
  
    vmd > atomselect macro oxygen
    name "O.*"
    
    vmd > atomselect macro hydrogen
    name "[0-9]?H.*"
    
    vmd > atomselect macro nitrogen
    name "N.*"
    
    vmd > atomselect macro carbon
    name "C.*" and not ion
  

  It was felt that, rather than slavishly follow VMD in this respect, Msys
  should try to get the correct answer.  Do you really want your "nitrogen"
  atom selection to include sodium (NA)?  If you do wish to recover the
  original VMD behavior, you can always redefine any of the atom selection
  macros as described below.

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

The Msys atom selection language
--------------------------------

For those not already familiar with the VMD atom selection language,
the following is a summary of the built-in atom selection facilities
of Msys.  The selection language can be extended on a per-System basis
in two ways, described in the subsequent sections.


Grammar
*******

The Msys atom selection grammar supports the following primitive types of
selections.

* Keyword selections: an attribute followed by one or more values, ranges,
  or regular expressions::

    name CA             # a single value
    resid 10 20 30      # multiple values
    resid 10 to 30      # a range of values
    name "C.*"          # regular expression, recognized by double quotes
    resid 5 8 to 10     # combination of single value and a range

* Singleword selections: A boolean attribute.  This includes both
  built-in singlewords such as ``protein``, as well as those defined
  as macros, such as ``hydrogen`` ::

    protein             # selects atoms identified as protein
    hydrogen            # defined by default as 'atomicnumber 1'

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
     
In comparisons, expressions can be formed in the following ways.

* Numeric literals, and keywords with numeric types::

    x < 3                       # atoms whose x coordinate is less than 3

* Composition of expressions using various mathematical operators, 
  parentheses, and functions::

    x + y * z < 3               # the usual precedence rules apply
    sqr(x)/36 + sqr(z)/125 < 1  # an ellipsoidal cylinder


Built-in selections
*******************

The following selection keywords are available:

  ================  =========== ===========================================
  keyword           type        definition
  ================  =========== ===========================================
  atomicnumber      integer     `Atom`.atomic_number
  element           string      Abbreviation for element `Atom`.atomic_number
  chain             string      `Chain`.name
  charge            float       `Atom`.charge
  fragment          integer     Connected residues will all have the same 
                                fragment id, except that the connection 
                                check will not follow disulfide bridges, 
                                identified as atoms whose name is "SG".
  index             integer     `Atom`.id
  mass              float       `Atom`.mass
  name              string      `Atom`.name
  numbonds          integer     `Atom`.nbonds
  resid             integer     `Residue`.resid
  residue           integer     `Residue`.id
  resname           string      `Residue`.name
  fragid            integer     `Atom`.fragid.  Connnected atoms will all 
                                have the same fragid.
  x, y, z           float       `Atom`.x, `Atom`.y, `Atom`.z, the position.
  vx, vy, vz        float       `Atom`.vx, `Atom`.vy, `Atom`.vz, the velocity.
  ================  =========== ===========================================


The following selection singlewords are available.  Note that these
are not implemented as macros and thus cannot be overridden or removed
by the user.

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


The following built-in macros are defined when a System is first created.
Users are free to override or delete them.

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
  bonded        numbonds > 0
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
  ion           resname AL BA CA Ca CAL CD CES CLA CL Cl CO CS CU Cu CU1 CUA HG IN IOD K MG MN3 MO3 MO4 MO5 MO6 NA Na NAW OC7 PB POT PT RB SOD TB TL WO4 YB ZN ZN1 ZN2
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

The atom selection language can be extended on a per-System basis with
macros.  A macro must be a single word, and cannot conflict with existing
selection keywords such as ``name``.  There are a number of pre-defined atom selection macros,
which you can list with `Selection.selection_macros`.  Other methods
in `System` let you view or change the definition of a macro, or remove
it altogether from the language.  Your changes to the selection macros
are saved in DMS files.

One use case for atom selection macros is when you have to work with
multiple related chemical systems with different atom selections for
corresponding functional groups. For example, the "active site" may
correspond to residues 32, 40, 48 for one chemical system, but residues
30, 31, 43, and 47 in another system.  If you define the atom selection
macro appropriately for each system and save it in the DMS file, you
will be able to simply select "active_site" when working with either
file and it will just work::

    mol.addSelectionMacro('active_site', 'chain A and resid 32 40 48')
    sel=mol.select('same residue as water and within 3 of active_site')

Atom selection macros can also be listed and updated using the 
``dms-macro`` command line tool.



