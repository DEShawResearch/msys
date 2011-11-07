
----------
Conversion
----------

mae2dms
-------
.. program:: mae2dms

.. describe:: mae2dms input.mae output.dms

   Converts an mae file to a dms file, preserving as much forcefield
   information as possible.

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


dms-dump
---------
.. program:: dms-dump

.. describe:: dms-dump file.dms [ options ]

   Writes a readable, grep-able summary of a dms file to stdout.

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

------------------
Structure building
------------------

dms-glue
--------
.. program:: dms-glue

.. describe:: dms-glue input.dms output.dms [-s selection]

   Finds a minimal set of "glue" bonds between the atoms in the selection,
   and writes those bonds to a "glue" table in the output file.

options:

.. cmdoption:: -s selection, --selection selection

   Selects atoms from the input dms file.  Default 'protein'.

.. cmdoption:: -v, --verbose

   Be chatty.

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


dms-solvate
-----------
.. program:: dms-solvate

.. describe:: dms-solvate input.dms output.dms [ options ]

   Adds a solvation box around the input structure.

dms-neutralize
--------------

.. program:: dms-neutralize

.. describe:: dms-neutralize input.dms output.dms [ options ]

   Replaces water molecules with ions in order to achieve a desired
   ion concentration.



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

