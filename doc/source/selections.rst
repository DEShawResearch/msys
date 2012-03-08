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

* ``atomicnumber`` (integer) `Atom`.atomic_number

* ``chain`` (string) `Chain`.name

* ``charge`` (float) `Atom`.charge

* ``fragment`` (integer) Connected residues will all have the same fragment
  id, except that the connection check will not follow disulfide bridges,
  identified as atoms whose name is "SG".

* ``index`` (integer) `Atom`.id

* ``mass`` (float) `Atom`.mass

* ``name`` (string) `Atom`.name

* ``numbonds`` (integer) `Atom`.nbonds

* ``resid`` (integer) `Residue`.resid

* ``residue`` (integer) `Residue`.id

* ``resname`` (string) `Residue`.name

* ``fragid`` (integer) `Atom`.fragid.  Connnected atoms will all have the same
  fragid.

* ``x``, ``y``, ``z`` (float) `Atom`.x, `Atom`.y, `Atom`.z, the position.

* ``vx``, ``vy``, ``vz`` (float) `Atom`.vx, `Atom`.vy, `Atom`.vz, the 
velocity.


The following selection singlewords are available.  Note that these
are not implemented as macros and thus cannot be overridden or removed
by the user.

* ``water`` -- atoms belonging to a residue containing the atomic number
  and bond structure of water, as well as those residues whose residue
  name is one of "H2O", "HH0", "OHH", "HOH", "OH2", "SOL", "WAT", "TIP", 
  "TIP2", "TIP3", "TIP4", or "SPC".

* ``hydrogen`` -- atomic number 1.

* ``backbone`` -- This singleword includes both protein backbone as well
  as nucleic backbone.  Protein backbone is identified by searching for
  atoms named "CA", "C", "O", and "N" in the same residue, and for atoms
  named "OT1", "OT2", "OXT", "O1", or "O2" that are bonded to one of the
  members of the first list.  If at least four such atoms are found, those
  atoms are identified as backbone.  Similarly, nucleic acid backbone atom
  names are P", "O1P", "O2P", "OP1", "OP2", "C3*", "C3'", "O3*", "O3'",
  "C4*", "C4'", "C5*", "C5'", "O5*", or "O5'"; or atoms named "H5T" or
  "H3T" bonded to a member of the first set.  At least four such atoms
  must be found in the same residue in order to be identified as backbone.

* ``protein`` -- residues containing protein backbone atoms.

* ``nucleic`` -- residues containing nucleic backbone atoms.


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


