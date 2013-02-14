
**************
DMS versioning
**************

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

---------------------
Schema changes in 1.7
---------------------

An **msys_ct** table was added, which stores arbitrary key-value pairs
for components of the full system.  The **particle** table has a new
column called *ct* which assigns the particle to a row of the **msys_ct**
table.  

Msys 1.6.x should be able to read msys 1.7.x files with no problems, since
the ct information will just be ignored with no effect on atom selections
or forcefield.

