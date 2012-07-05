
----------
Conversion
----------

mae2dms
-------
.. program:: mae2dms

.. describe:: mae2dms input.mae output.dms

   Converts an mae file to a dms file, preserving as much forcefield
   information as possible.

options:

.. cmdoption:: --ignore-unrecognized   

   skip unrecognized ffio_ff subblocks

*mae2dms* converts mae files to dms files.  Atom order and forcefield
information are all preserved, though there are likely to be differences
in the precise values of the forces computed from the file due differences
in how force terms are represented.

The *m_chain_name* field is used to name the chains created by Msys.  If
the *m_pdb_segment_name* field is present, it will be mapped to the
custom particle property *segid*.  Thus, atoms in the same residue may
have different segids, but by construction will always have the same chain.

dms2mae
-------
.. program:: dms2mae

.. describe:: dms2mae input.dms output.mae

    Converts a dms file to an mae file, preserving as much forcefield
    information as possible.

*dms2mae* is the inverse of *mae2dms*, though with somewhat more restricted
range of supported forcefield terms. 

-----------
Information
-----------

dms-info
--------
.. program:: dms-info

.. describe:: dms-info [ options] [ dms files ]

   Writes a summary of the atom and forcefield information of a dms file.

*dms-info* provides a summary of the contents of a dms file.  Its output
includes:

 * total number of atoms, bonds, residues, and chains

 * global cell size

 * force tables organized by category (bond, constraint, virtual, etc).

 * number of atoms each that can be selected as protein, lipid, ions, water,
   or none of the above.

::

  > dms-info leuTaa_leu_POPC.dms 
  ---------------------------------------------------------------------------
  leuTaa_leu_POPC.dms
  
  Structure:
         Atoms:     8982
         Bonds:     8948
      Residues:      675
        Chains:        1
   Global cell: (81.659309, 0.0, 0.0)
                (0.0, 81.659309, 0.0)
                (0.0, 0.0, 82.648071)
  
       Protein:     8225 atoms,      512 residues,        1 chains
         Lipid:      268 atoms,        2 residues,        1 chains
          Ions:        0 atoms,        0 residues,        0 chains
         Water:      477 atoms,      159 residues,        1 chains
         Other:       12 atoms,        2 residues,        1 chains
  
  Bond Tables:
         alchemical_angle_harm: 3 sites,      2 params,      2 terms
      alchemical_dihedral_trig: 4 sites,      0 params,      0 terms
      alchemical_improper_harm: 4 sites,      3 params,      3 terms
       alchemical_pair_12_6_es: 2 sites,     51 params,     51 terms
       alchemical_stretch_harm: 2 sites,      1 params,      1 terms
                    angle_harm: 3 sites,    130 params,  15664 terms
                 dihedral_trig: 4 sites,     71 params,  22942 terms
                 improper_harm: 4 sites,     16 params,   1200 terms
                  pair_12_6_es: 2 sites,  22158 params,  22412 terms
                  stretch_harm: 2 sites,     89 params,  12235 terms
           torsiontorsion_cmap: 8 sites,      5 params,    511 terms
  
  Constraint Tables:
                constraint_ah1: 2 sites,     10 params,   1652 terms
                constraint_ah2: 3 sites,      5 params,    703 terms
                constraint_ah3: 4 sites,      4 params,    424 terms
                constraint_hoh: 3 sites,      1 params,    159 terms
  
  Exclusion Tables:
                     exclusion: 2 sites,      0 params,  47234 terms
  
  Nonbonded Tables:
          alchemical_nonbonded: 1 sites,     50 params,     22 terms
                     nonbonded: 1 sites,     50 params,   8982 terms
  
  Nonbonded Info:
             vdw_funct: vdw_12_6
              vdw_rule: arithmetic/geometric
  
          Glue:
  
  Auxiliary Tables:
                         cmap1: 3 properties,    576 rows
                         cmap2: 3 properties,    576 rows
                         cmap3: 3 properties,    576 rows
                         cmap4: 3 properties,    576 rows
                         cmap5: 3 properties,    576 rows
                   dms_version: 2 properties,      1 rows
                    forcefield: 3 properties,      3 rows
                   viparr_info: 3 properties,      1 rows
  
  Provenance:
  

dms-dump
---------
.. program:: dms-dump

.. describe:: dms-dump file.dms [ options ]

   Writes a readable, line-based (i.e., grep-able) summary of a dms file
   to stdout.

options:

.. cmdoption:: --without-provenance

   Don't print the provenance section of the dms file.

.. cmdoption:: --without-groups

   Don't print columns in the particle table beginning with ``grp_``.

.. cmdoption:: --without-forcefield

   Don't print the forcefield information section of the dms file.


*dms-dump* generates a textual representation of a dms file that can be
understood by humans and compared to the output of another dms file.
A certain amount of canonicalization is applied to the contents of the dms
file in order to make this happen:

 * all floating point values are rounded to six decimals;

 * all force terms are sorted by particle id, i.e. p0, p1, ...

 * all columns are sorted alphabetically;

 * rather than printing the nbtype of each particle, the nonbonded parameters
   of each particle are dumped with a "nonbonded" label at the start of
   each line.


dms-diff
--------
.. program:: dms-diff

.. describe:: dms-diff file1.dms file2.dms

  Writes a Unix diff of the dms files ``file1.dms`` and ``file2.dms`` to
  standard output.  The environment variable ``DMSDIFF`` can be used to
  specify an alternate file comparison utility.


------------------
Basic Manipulation
------------------

dms-frame
---------
.. program:: dms-frame

.. describe:: dms-frame input.dms output.dms [ options ]

   Extract a frame from a trajectory into a dms file, and perform optional
   centering and periodic wrapping.

options:

.. cmdoption:: -i input, --input-path input

   Input coordinate/trajectory file.

.. cmdoption:: --input-type type

   File type for input file; default 'dtr'.

.. cmdoption:: -t time, --time time

   Selected frame time for input coordinates

.. cmdoption:: -n index, --index index

   Selected frame index for input coordinates

.. cmdoption:: --zero-velocities

   Use zero velocities instead of reading from frame.

.. cmdoption:: -g glue, --glue glue

   Glue atom selections (can be specified multiple times)

.. cmdoption:: -c centersel, --center centersel

   Center atoms in the given selection

.. cmdoption:: --wrap

   Apply periodic wrapping (implied by --center)


*dms-frame* reads coordinates from a coordinate or trajectory file and
copies them into a DMS file.  Periodic wrapping may also be applied
to the coordinates.  If the input coordinate file is a trajectory containing
multiple frames, the frame may be selected either with a time, in which case
the frame whose time is closest to the provided value will be selected,
or with a frame index, which can be negative in the usual Python sense
(i.e. index -1 chooses the last frame).  Either time or index may be 
provided, but not both.

If a centering selection is provided, the center will be calculated from
the input coordinates after applying any glue or periodic wrapping.  The
selection will be centered on the origin, and the rest of the system will
be wrapped so as to fit in the periodic cell.

Velocities will be copied from the input frame if they are present; if
not, the velocities in the input dms file will be used.  Specifying 
`--zero-velocities` makes the velocities zero in the output file.


dms-select  
----------
.. program:: dms-select

.. describe:: dms-select input.dms [ options ]

   Write or append a selection from ``input.dms`` to an output dms file.

options:

.. cmdoption:: -s selection, --selection selection

   Select atoms from the input dms file.

.. cmdoption:: -o output.dms, --output output.dms

   Write the selected atoms to ``output.dms``.

.. cmdoption:: -a output.dms, --append output.dms

   Append the selected atoms to ``output.dms``.

.. cmdoption:: -v, --verbose

   Print information about the selected atoms to stdout.

*dms-select* takes the selected atoms in ``input.dms`` and either writes
or appends them to ``output.dms``, depending on the supplied options.


dms-set
-------
.. program:: dms-set

.. describe:: dms-set input.dms output.dms [ options ] [ updates ]

   Updates atom, residue, chain, and/or table properties of the particles
   in input.dms; writes result to output.dms.


options:

.. cmdoption:: -s selection, --selection selection

   Selects atoms from the input dms file to update.

update format:

.. cmdoption:: atomprop=FOO

   Change the value of atom property ``atomprop`` to FOO.

.. cmdoption:: residue.resprop=BAR

   Change the value of residue property ``resprop`` to BAR.
  
.. cmdoption:: chain.chainprop=BAZ

   Change the value of chain property ``chainprop`` to BAZ.

.. cmdoption:: table.tableprop=XYZ

   Change the property `'tableprop`` in table ``table`` to XYZ.


*dms-set* creates a new dms file with modifications to the atom, residue,
chain, or table properties.  Multiple updates may be specified, in which
case they will be applied in the order they are given.  If an atom selection
is provided, it is evaluated before any of the updates are applied.

Updates to residues and chains are applied to every residue or chain
with `at least one atom` in the selection.  Updates to tables are applied
to terms whose atoms are `completely contained` in the selection.  Since
the update is specified in terms of an atom selection, the order of atoms
in the terms is irrelevant in determining whether a term is affected by
a update.

Example: Change the name CD1 atoms in LYS residues to CD.

   ``dms-set input.dms output.dms -s "resname LYS and name CD1" name=CD``


Example: Change the stretch term force constant to 0 for a pair of atoms
with ids 32 and 42.  As described above, this would not affect the stretch
terms involving atoms 32 or 42 with any other atom; only the term involving
both atoms. 

   ``dms-set input.dms output.dms -s "index 32 42" stretch_harm.fc=0.0``


dms-macro
---------

.. program:: dms-macro

.. describe:: dms-macro system.dms [ options ]

   List and modify the set of atom selection macros in a dms file.

options:

.. cmdoption:: -d macro, --delete macro

   Remove `macro` from the system's selection macro list.

.. cmdoption:: -m macro="ATOM SELECTION", --macro macro="ATOM SELECTION"

   (Re)define macro named `macro` to the given "ATOM SELECTION".

.. cmdoption:: -l

   Print the macros and their defintions in the system.

.. cmdoption:: -o output.dms

   Write out the system with modifications to the macros to output.dms

*dms-macro* is used to query and modify the set of atom selection
macros available in the given dms file.  More information about macros
may be found the atom selections section of the documentation.


------------------
Structure building
------------------

dms-glue
--------
.. program:: dms-glue

.. describe:: dms-glue input.dms output.dms [-s selection] [--replace]

   Finds a minimal set of "glue" bonds between the atoms in the selection,

options:

.. cmdoption:: -s selection, --selection selection

   Selects atoms from the input dms file.  Default 'protein'.

.. cmdoption:: --replace

   Replace any existing glue pairs in the system.

.. cmdoption:: -v, --verbose

   Be chatty.


*dms-glue* finds pairs of atoms which connect the disjoint sets of atoms
in a selection together.  The default algorithm uses a heuristic to find
nearby atoms in different connected subsets to serve as the glue pairs.


dms-grease
----------
.. program:: dms-grease
  
.. describe:: dms-grease input.dms lipid.dms output.dms [ options ]

   Adds a lipid bilayer around a solute.

.. cmdoption:: --structure-only

   Load only the structure part of input.dms and lipid.dms, not the forcefield

.. cmdoption:: -t thickness, --thickness thickness

   Minimum distance from outer edge of membrane to input structure

.. cmdoption:: -x xsize, --xsize xsize

   Size of membrane along x dimension.  Overrides --thickness.

.. cmdoption:: -y ysize, --ysize ysize

   Size of membrane along y dimension.  Overrides --thickness.

.. cmdoption:: -c chain, --chain chain

   Chain name of constructed bilayer

.. cmdoption:: --square

   Ensure xsize and ysize are equal to max(xsize, ysize)

.. cmdoption:: -v, --verbose

   Be chatty.

dms-grease builds a new chemical system consisting of the input system
plus a lipid bilayer constructed by tiling *lipid.dms* in the x-y plane.
If the *input.dms* is given as "-", then a pure membrane will be built.

An error will be encountered if only one of *input.dms* and *lipid.dms* 
have forcefield information; this is because Msys refuses to write DMS
files that have only partial information for the nonbonded atom types.
If you don't have forcefield information for one of the input files,
use the *--structure-only* option to ignore the forcefield information
in the one that does.

The global cell of the new system will be orthorhombic and have x and
y dimensions given by the specified size of the membrane, and z dimension
given by the input structure or the lipid membrane template, whichever is
greater.


dms-thermalize
--------------

.. program:: dms-thermalize

.. describe:: dms-thermalize input.dms output.dms [ options ]

   Assign Boltzmann-sampled velocities to the atoms.  Atoms with zero mass
   will get zero velocity.

.. cmdoption:: -t TEMPERATURE, --temperature TEMPERATURE

   Sample Boltzmann distribute with given temperature in Kelvin.

.. cmdoption:: -s SEED, --seed SEED

   Use the given random seed, default 1, or 'random' to get a random random 
   seed.


dms-posre
---------

.. program:: dms-posre

.. describe:: dms-posre input.dms output.dms [ options ]

   Assign harmonic position restraints to selected atoms.  


.. cmdoption:: -f FORCE_CONSTANT

   force constant in PEAK units (kcal/mol/A^2)

.. cmdoption:: -x FORCE_CONSTANT

   force constant along x axis in PEAK units (kcal/mol/A^2)

.. cmdoption:: -y FORCE_CONSTANT

   force constant along y axis in PEAK units (kcal/mol/A^2)

.. cmdoption:: -z FORCE_CONSTANT

   force constant along z axis in PEAK units (kcal/mol/A^2)

.. cmdoption:: -s selection, --selection=selection

   Add/replace position restraint for selected atoms

.. cmdoption:: --replace

   Remove all existing position restraints.

.. cmdoption:: --quiet

   Turn off chattiness


`dms-posre` adds position restraints to a dms file, using the existing atom
positions for the reference positions of the restraints.  If ``--replace``
is specified on the command line, any existing restraints will be replaced
by the new set.  Otherwise, atoms that are already restrained in the existing
file will be restrained using the newly provided force constraints::

  # Add position restraints to backbone atoms with a force constant of 0.2
  dms-posre input.dms out1.dms -s "backbone" -f 0.2

  # Restrain CA atoms with a force constant of 0.3
  dms-posre out1.dms out2.dms -s "name CA" -f 0.3

dms-override-vdw
----------------
.. program:: dms-override-vdw

.. describe:: dms-override-vdw input.dms output.dms [ options ]

   Override vdw interactions between selected atoms.

options:

.. cmdoption:: --sigma sigma

   Vdw sigma

.. cmdoption:: --epsilon epsilon

   Vdw epsilon

.. cmdoption:: --selection0 selection
 
   Atom selection for the first group

.. cmdoption:: --selection1 selection

   Atom selection for the second group

*dms-override-vdw* changes the vdw interaction between two specified groups
of atoms to the specified values of sigma and epsilon.  All options (sigma,
epsilon, selection0, selection1) are required, and the selection groups must
not be empty.  

Currently, the vdw functional form of the DMS file must be "vdw_12_6".  

This tool uses the `nonbonded_combined_param` table in the DMS file to store
the overrides and therefore should not be used with versions of Anton
software older than 2.9.2  

dms-scale-vdw
-------------
.. program:: dms-scale-vdw

.. describe:: dms-scale-vdw input.dms output.dms [ options ]

   Scale vdw interactions between selected atoms.

options:

.. cmdoption:: -s scale_sigma, --scale-sigma scale_sigma

   scale factor for sigma, default 1.0

.. cmdoption:: -e scale_epsilon, --scale-epsilon scale_epsilon

   scale factor for epsilon, default 1.0

.. cmdoption:: -l selection, --ligand selection

   atom selection groups (specify multiple)

*dms-scale-vdw* scales the vdw interactions between multiple groups of atoms.
The vdw interactions between each ligand group will be scaled by the
specified amount.  As many ligands may be specified as desired, though
different implementations on Desmond and Anton may in practice limit the
number possible

Currently, the vdw functional form of the DMS file must be "vdw_12_6".  

This tool uses the `nonbonded_combined_param` table in the DMS file to store
the overrides and therefore should not be used with versions of Anton
software older than 2.9.2  


------------------------
Free Energy Perturbation
------------------------

dms-uncharge
------------
.. program:: dms-uncharge

.. describe:: dms-uncharge input.dms output.dms [ options ]

   Create an alchemical dms file with selected atoms uncharged in the B state.

.. cmdoption:: -s selection, --selection selection

   Uncharge only atoms in selection


dms-alchemical
--------------
.. program:: dms-alchemical

.. describe:: dms-alchemical input.dms output.dms atom.map C.dms 

   Create an alchemical system from A and B states and a map between them.


The *atom.map* file should consist of lines with two 1-based indices,
the first referring to atoms in the A state and the second to atoms in
the B state.  Either the A or B index may be negative, indicating that
the corresponding atom has no analog in the other state.  The mapping
must reference the first Na atoms in the A state and Nb atoms in the B
state, where Na need not equal Nb.  

The generated alchemical system C will have N alchemical atoms, where N
is the number of lines in *atom.map*.   Atoms and force terms in the A state
not referenced by the atom map will be appended to the structure; unreferenced
atoms and force terms in the B state will be ignored.

----------
Validation
----------


dms-validate
------------
.. program:: dms-validate

.. describe:: dms-validate input.dms [ options ]

    Perform various sanity checks on a chemical system.

.. cmdoption:: --strict

    Also perform strict checks.

.. cmdoption:: --desmond

    Also perform Desmond-specific checks.

.. cmdoption:: --verbose

    Be verbose.

`dms-validate` flags conditions that are likely to be errors in a chemical
system.  The set of "basic" checks are always performed; additional checks
can be enabled using various command line flags. 

The set of basic checks comprise the following:

 * nonbonded: if a nonbonded table exists, every particle must have a 
   nonbonded param assignment.


The set of strict checks comprise the following items.  Note that it
is certainly possible for a valid simulation to be performed using a
system that passes none of its strict checks!  However, it may be worth
investigating why a system fails theses checks.

 * constraints: the system must have constraint terms.  

 * consistent masses: Particles with equal atomic number must have equal mass.
   Pseudo particles (those with atomic number zero) are excluded from the
   check.

 * sparsify: every 1-4 bond (i.e., pair of atoms separated by three 
   distinct bonds) must be included in the exclusion table.

Desmond-specific checks:

 * bonded terms: check that neither the exclusion table nor any table
   in the bond_term metable contains terms whose atoms are not connected
   through the bond table.  


