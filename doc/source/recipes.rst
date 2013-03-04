

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


