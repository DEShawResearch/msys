
*************************
Notes on MAE file support
*************************

As of this writing (January 2013), the MAE file format still sees 
common use within DESRES and beyond.  The following notes are intended
to guide the expectations of users of MAE files when interacting
with msys and msys-derived tools such as viparr.  In particular,
we take note of MAE file contents that cannot be represented
within Msys at all, due to specific aspects of the design of Msys
that conflict with MAE files, and vice versa: not all Msys systems
can be readily serialized to an MAE file.

Multiple ct blocks
------------------

When there is more than one "f_m_ct" block in a MAE file, it can mean
either that the file contains multiple structures, or that the individual
blocks are separated for bookkeeping purposes but should be considered
to be a single physical system.  The Desmond and Anton software takes
the latter point of view, and it is therefore necessary to merge or
otherwise reconcile redundant information in the ct blocks when treating
them as a single system.

The LoadMany function in Python (LoadIterator in C++) loads each ct block
separately into a new System.

The LoadMAE function in Python (ImportMAE in C++) loads all ct blocks in
an MAE file and merges them into a single System.  The following 
requirements must be satisfied on the input blocks:

 * The global cell parameters (r_chorus_box_{abc}{xyz}) of the _last_
   ct block takes precedence.

 * The nonbonded functional form in ffio_vdwtypes, and the combining
   rule, must be the same in all blocks.  If the rule is left blank,
   rules in prior blocks takes precedence.


Ct-level attributes
-------------------

When using LoadMany, the "name" attribute of each System will be read from
the "m_title" field of the corresponding ct block.  LoadMAE will set the
System name to the path name of the MAE file.

Msys reads and writes a ct-level array called "msys_provenance" in order
to track the history of its files.

No other Ct-level attributes are stored or read by Msys, other than 
global cell data.

Pair terms
----------


Virtual sites
-------------


VDW overrides (NBFix)
---------------------


Multiple nonbonded tables
-------------------------



