################
The Msys library
################

Msys is a set of libraries and tools for:

* reading and writing chemical structure and trajectory files;
* editing and combining structures and forcefields.

It comes with command line tools which automate some common tasks, as well as C++ and Python
interfaces which offer more power and flexibility.

To illustrate what msys can do, use the dms-info tool to inspect the content of a chemical
system. Here's the result of running dms-info on a typical pdb file::

    $ dms-info 5ak4.pdb
    ---------------------------------------------------------------------------
    5ak4.pdb
    Structure :
           Atoms:     4978 (0 pseudo)
           Bonds:     4592
        Residues:     1129
          Chains:        2
     Global cell:     92.339        0.0        0.0
                  -46.16949999999999 79.96791976005129        0.0
                         0.0        0.0    245.024
    
         Protein:     4343 atoms,      546 residues,        1 chains
           Lipid:        0 atoms,        0 residues,        0 chains
         Nucleic:        0 atoms,        0 residues,        0 chains
            Ions:        0 atoms,        0 residues,        0 chains
           Water:      578 atoms,      577 residues,        1 chains
           Other:       57 atoms,        6 residues,        2 chains
    
    Ct composition:
           1)                              4978 atoms
    
    Nonbonded Info:
               vdw_funct:
                vdw_rule:
    
    Auxiliary Tables:
    
    Forcefields:
    
    Provenance:


Msys organizes the atoms in the system into residues, chains, and components. It also reports
commonly used types of chemical matter, including protein and water. The bonds reported by
dms-info for this pdb file were inferred by msys since pdb files typically lack connectivity
information for most of the structure.

One of the most useful features of msys is its atom selection language. The dms-select tool
illustrates how you can use atom selections to pull out subsets of the system::

    $ dms-select 5ak4.pdb -s "same fragid as protein" -o 5ak4-protein.pdb
    Loading 5ak4.pdb
    Cloning "same fragid as protein"
    Selected 2 chain(s), 548 residue(s), 4346 atom(s)
    Saving 5ak4-protein.pdb
    $ dms-info 5ak4-protein.pdb
    ---------------------------------------------------------------------------
    5ak4-protein.pdb
     
    Structure :
           Atoms:     4346 (0 pseudo)
           Bonds:     4505
        Residues:      548
          Chains:        1
     
     Global cell:     92.339        0.0        0.0
                  -46.16949999999999 79.96791976005129        0.0
                         0.0        0.0    245.024
     
         Protein:     4343 atoms,      546 residues,        1 chains
           Lipid:        0 atoms,        0 residues,        0 chains
         Nucleic:        0 atoms,        0 residues,        0 chains
            Ions:        0 atoms,        0 residues,        0 chains
           Water:        2 atoms,        1 residues,        1 chains
           Other:        1 atoms,        1 residues,        1 chains
     
    Ct composition:
           1)                              4346 atoms
     
    Nonbonded Info:
               vdw_funct:
                vdw_rule:
     
    Auxiliary Tables:
     
    Forcefields:
     
    Provenance:

The atom selection "same fragid as protein" includes protein atoms as well as anything connected
to protein via bonds.

Msys can read and write coordinate sets (frames) from a trajectory, process them in various
ways, and write them into new chemical system files. For example, here is how you could use
the dms-frame tool to copy the last frame of a trajectory into a chemical system, preserving
the existing atom data and applying a centering operation on the coordinates:

.. code-block::

    $ dms-frame system.dms system-relax.dms -i relax.dtr -n -1 -c protein

Here, `system.dms` and `system-relax.dms` are the input and output chemical system files, and
`relax.dtr` holds trajectory data.


Reading and writing chemical systems from Python is easy::

  import msys

  # Load the entire contents of a DMS file
  dms=msys.LoadDMS('system.dms')

  # Import an MAE file, performing conversions on its forcefield data
  mae=msys.LoadMAE('system.mae')

You can also create a new `System` from scratch::

  # mol = msys.CreateSystem()

A `System` resides entirely in memory; changes to the `System` will not
persist until/unless you write it back out to a file::

  # Save the system as a DMS file
  msys.SaveDMS(dms, 'output.dms')

  # Export to MAE file
  msys.SaveMAE(dms, 'output.mae')

The rest of this guide goes into more detail on the atom selection language, how msys represents
forcefields, tools for working with trajectories and coordinates, and more.

.. toctree::
   :maxdepth: 2

   overview
   selections
   forcefield
   wrapping
   python
   tools
   recipes
   dms
   release_notes


