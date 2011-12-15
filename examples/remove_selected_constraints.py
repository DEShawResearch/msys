#!/usr/bin/env python2.7

'''
Remove all constraints for atoms in a particular atom selection, and also
update the 'constrained' property of the stretch_harm terms for the
stretch terms that are no longer constrained.
'''

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

