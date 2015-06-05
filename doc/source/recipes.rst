
Adding artificial bonds
-----------------------

Msys can add force terms to a system::

    import msys, sys
    ifile, ofile = sys.argv[1:]
    mol=msys.Load(ifile)
    T=mol.table('stretch_harm')
    P=T.params
    param=P.addParam()
    param['fc']=32
    param['r0']=1.0
    T.addTerm([mol.atom(0), mol.atom(1)], param)
    # ...
    msys.Save(mol, ofile)



Adding energy groups
--------------------

Desmond and Anton use the "energy_groups" atom property to assign atoms to
energy groups::

  mol = msys.Load('system.dms')

  # add an integer property.  The default value is zero.  It's a no-op
  # if the property already exists, and an error if it exists but has a
  # different type.
  mol.addAtomProp('grp_energy', int)        

  # assign protein to energy group 1
  for a in mol.select('protein'):
    a['grp_energy'] = 1

  # save the result
  msys.SaveDMS(mol, 'system_engrp.dms')


Remove selected constraints
---------------------------

Remove all constraints for atoms in a particular atom selection, and also
update the 'constrained' property of the stretch_harm terms for the
stretch terms that are no longer constrained::

    import sys, msys
    ifile, ofile=sys.argv[1:]
    seltext='hydrogen and withinbonds 1 of (backbone and nitrogen)'
    
    print "reading from %s" % ifile
    mol=msys.LoadDMS(ifile)
    atoms=set(mol.select(seltext))
    residues=set(a.residue for a in atoms)
    
    print "found %d atoms in selection from %d residues" % (len(atoms), len(residues))
    
    removed_constraints=list()
    
    for table in mol.tables:
        if table.category != 'constraint': continue
        print "doing table", table.name
        for term in table.terms:
            for a in term.atoms:
                if a in atoms:
                    removed_constraints.append(set(term.atoms))
                    term.remove()
                    break
    
    print "removed %d constraints" % len(removed_constraints)
    
    stretch=mol.table('stretch_harm')
    if 'constrained' in stretch.term_props:
        print "finding constrained stretch_harm terms"
    
        uncons=0
        for term in mol.table('stretch_harm').terms:
            if not term['constrained']: continue
            s=set(term.atoms)
            # one of the atoms in the term must be in the original selection
            if not s.intersection(atoms):
                continue
    
            # mark as unconstrained if all atoms in term overlap a constraint
            for cons in removed_constraints:
                if not s.difference(cons):
                    term['constrained']=0
                    uncons += 1
                    break
    
    print "unconstrained %d stretch_harm terms" % uncons
    
    print "Saving to %s" % ofile
    msys.SaveDMS(mol, ofile)


Canonicalize position restraint terms
-------------------------------------

Some posre_harm tables have x0, y0, z0 as param properties rather than
term properties.  This script enforces the convention that these properties
ought to be term properties::

    def canonicalize(mol):
        posre=mol.table('posre_harm')
        props=set(('x0', 'y0', 'z0'))
        if props.issubset(set(posre.term_props)):
            print "Already canonical!"
            return
    
        if not props.issubset(set(posre.params.props)):
            print "Missing %s from posre params!" % (props,)
            exit(1)
    
        print "File is not canonical!  Fixing..."
        posre.name = '__posre_harm_old__'
        newposre=mol.addTableFromSchema('posre_harm')
        for t in posre.terms:
            p = newposre.params.addParam()
            p['fcx'] = t['fcx']
            p['fcy'] = t['fcy']
            p['fcz'] = t['fcz']
            t2 = newposre.addTerm(t.atoms, p)
            t2['x0'] = t['x0']
            t2['y0'] = t['y0']
            t2['z0'] = t['z0']
        posre.remove()
        newposre.coalesce()

    def main():
        import sys
        ifile, ofile = sys.argv[1:]
        mol=msys.LoadDMS(ifile)
        canonicalize(mol)
        mol = mol.clone()
        msys.SaveDMS(mol, ofile)
        
    if __name__=="__main__": main()

Processing multi-entry files (e.g. SDF files)
---------------------------------------------

To iterate over each structure in an SDF file, use msys.LoadMany.
The LoadMany function is a generator, so you should iterate over its
results rather than simply calling it.

Each result returned by LoadMany is a System with one ct.  You'll
need to access the ct member of the System in order to view or modify
the data values associated with each entry.  

To write entries back out to a new SDF file after filtering or modifying
them, use msys.SaveSDF.  It's most efficient to create your own file
object in Python, and write the string returned by SaveSDF to that file.

When entries in the SDF file cannot be parsed, msys skips the next entry,
and msys.LoadMany returns None for the offending entry.  You should check
for None in the return values of msys.LoadMany and skip them if that makes
sense for your script.

Here is an example snippet which reads each entry, filters by atom count,
modifies a data property, removes another data property, and writes the
results to another file::

    def process_sdf(ifile, ofile):
        fp = open(ofile, 'w')
        for i, mol in enumerate(msys.LoadMany(ifile)):
            # skip entries which couldn't be parsed
            if mol is None:
                print "Warning, skipping entry %d" % (i+1)
                continue
            # filter systems with fewer than 5 atoms
            if mol.natoms < 5:
                continue
            ct = mol.ct(0)
            # update 'THE_SCORE' property.  Note that vlaues returned by
            # get may be float, int, or string.
            score = ct.get('THE_SCORE', 0.0)
            score += 5.0
            ct['THE_SCORE'] = score
            # remove 'USELESS' property
            if ct.get('USELESS') is not None:
                del ct['USELESS']
            # write the entry back out
            fp.write(msys.SaveSDF(mol, None)


Processing large SDF files
--------------------------

If you have large sdf files with many thousands of entries, you may
benefit from using a set of functions specialized for SDF files.
The new functions are around 10x faster at reading SDF files and
20x faster for writing.  However, there is no facility for modifying
the molecular structures of each entry, though you can inspect and
modify the data values.  Also, the data values are always returned as
strings, so you must case them to appropriate types if you wish to
manipulate them as integers or floats.

The new functions are named ScanSDF and FormatSDF.  Here a snippet
which performs the same actions as the process_sdf function in the
previous example, using the new functions::

    def process_sdf_fast(ifile, ofile):
        fp = open(ofile, 'w')
        for i, mol in enumerate(msys.ScanSDF(ifile)):
            # skip entries which couldn't be parsed
            if mol is None:
                print "Warning, skipping entry %d" % (i+1)
                continue
            # filter systems with fewer than 5 atoms
            if mol.natoms < 5:
                continue
            # update 'THE_SCORE', coverting the existing value to float
            score = mol.get('THE_SCORE', 0.0)
            ct['THE_SCORE'] = float(score) + 5.0
            # remove 'USELESS' property
            if ct.get('USELESS') is not None:
                del ct['USELESS']
            # write the entry back out
            fp.write(msys.FormatSDF(mol)
        
