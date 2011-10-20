
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

**Important note on MAE import:**  Msys uses a single
identifier called *chain* to represent groups of residues.  However,
MAE files, as well as PDB files, often contain both chain and segment
identifiers (segid).  an importer for Msys must somehow infer which one
to use for a given atom.  The approach taken by all the importers in
Msys is to use the *segid* whenever it contains non-whitespace,
and otherwise fall back to the *chain*.  The rationale is two-fold:
(1) *segid* is less commonly used than *chain*, so if it is present in the
file, it's likely to be what the user intended; (2) *segid* is allocated
more characters by the PDB format as well as programs such as VMD, so if
a file for some reason needed to include both segid and chain, the
importers here use *segid* so that longer and potentially more descriptive
names can be used.

dms2mae
-------
.. program:: dms2mae

.. describe:: dms2mae input.dms output.mae

    Converts a dms file to an mae file, preserving as much forcefield
    information as possible.

*dms2mae* is the inverse of *mae2dms*, though with somewhat more restricted
range of supported forcefield terms. 


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


