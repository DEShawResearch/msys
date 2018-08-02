'''
Add position restraints to a dms file, using the existing atom
positions for the reference positions of the restraints.  If ``--replace``
is specified on the command line, any existing restraints will be replaced
by the new set.  Otherwise, atoms that are already restrained in the existing
file will be restrained using the newly provided force constraints::

  # Add position restraints to backbone atoms with a force constant of 0.2
  dms-posre input.dms out1.dms -s "backbone" -f 0.2

  # Restrain CA atoms with a force constant of 0.3
  dms-posre out1.dms out2.dms -s "name CA" -f 0.3

  # Remove all position restraints:
  dms-posre input.dms output.dms --replace
  ## or:
  dms-posre input.dms output.dms -s none --replace
'''
from __future__ import print_function

import msys
import math

def apply(mol, atoms, fcx, fcy, fcz, replace=False):
    ''' add position restraints to atoms '''

    table=mol.addTableFromSchema('posre_harm')
    index={}
    if replace:
        for t in table.terms: t.remove()
    else:
        for t in table.terms:
            index[t.atoms[0]] = t

    param=table.params.addParam()
    param['fcx']=fcx
    param['fcy']=fcy
    param['fcz']=fcz

    for a in atoms:
        t=index.get(a)
        if t is None:
            t=table.addTerm([a])
        t.param = param
        t['x0']=a.x
        t['y0']=a.y
        t['z0']=a.z

    return table.nterms

def main():
    import optparse
    parser = optparse.OptionParser(__doc__)
    parser.add_option('-f', type='float', default=0.0,
            help='force constant in PEAK units')
    parser.add_option('-x', default=None,
            help='force constant along x axis in PEAK units')
    parser.add_option('-y', default=None,
            help='force constant along x axis in PEAK units')
    parser.add_option('-z', default=None,
            help='force constant along x axis in PEAK units')
    parser.add_option('-s', '--selection', default='none',
            help="selection for restrained atoms")
    parser.add_option('--replace', action="store_true", default=False,
            help='replace any existing harmonic position restraints')
    parser.add_option('--reference-structure', default=None,
            help='get equilibrium positions from reference structure')
    parser.add_option('--reference-selection', default=None,
            help='selection for reference structure [default: selection]')
    parser.add_option('--max-distance', default=None,
            help='abort if max restraint distance exceeds value in Angstroms')
    parser.add_option('--max-energy', default=None,
            help='abort if max restraint energy exceeds value in kcal/mol')

    parser.add_option('-q', '--quiet', action='store_true', default=False,
            help="Don't be chatty")

    opts, args = parser.parse_args()
    if len(args)!=2:
        parser.error("incorrect number of arguments")

    if not opts.quiet: print("Loading input file <%s>" % args[0])
    mol=msys.LoadDMS(args[0])

    atoms=mol.select(opts.selection)
    if not opts.quiet: print("Adding restraints to %d atoms" % len(atoms))

    fc=float(opts.f)
    fcx = fc if opts.x is None else float(opts.x)
    fcy = fc if opts.y is None else float(opts.y)
    fcz = fc if opts.z is None else float(opts.z)
    if not opts.quiet: print("Using force constant (%s, %s, %s)" % (
            fcx,fcy,fcz))

    ref=opts.reference_structure
    if ref:
        if not opts.quiet: print("Loading reference file <%s>" % ref)
        ref_mol = msys.LoadDMS(ref)
        if not opts.quiet: print("Restraining atoms to reference positions")
        if opts.reference_selection:
            ref_selection = opts.reference_selection
        else:
            ref_selection = opts.selection
        ref_sel = ref_mol.clone(ref_selection)
        cur_sel = mol.clone(opts.selection)
        assert len(ref_sel.atoms) == len(atoms)
        for a, p in zip(atoms, ref_sel.positions):
            a.pos = p

        displacements = ref_sel.positions - cur_sel.positions
        E_restr = ((displacements**2)*(fcx, fcy, fcz)).sum(-1)
        index = E_restr.argmax()
        dist = math.sqrt(sum(displacements[index]**2))
        a = atoms[index]
        r = a.residue
        if not opts.quiet:
            a_name = 'chain %s %s%d:%s' % (r.chain.name,r.name,r.resid,a.name)
            fmt = "Most energetic restraint: %s (%.1f kcal/mol, %.1f Angstroms)"
            print(fmt % (a_name, E_restr.max(), dist))
        if opts.max_distance is not None:
            if dist > float(opts.max_distance):
                print("ERROR: Maximum distance exceeded.", file=sys.stderr)
                return 1
        if opts.max_energy is not None:
            if E_restr.max() > float(opts.max_energy):
                print("ERROR: Maximum energy exceeded.", file=sys.stderr)
                return 1

    else:
        if not opts.quiet: print("Restraining atoms to current positions")
    n=posre.apply(mol, atoms, fcx, fcy, fcz, replace=opts.replace)
    if ref:
        for a, p in zip(atoms, cur_sel.positions):
            a.pos = p
    if not opts.quiet: print("Now have restraints on %d atoms" % n)

    if not opts.quiet: print("Writing DMS file <%s>" % args[1])
    msys.SaveDMS(mol,args[1])

