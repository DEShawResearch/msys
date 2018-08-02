from __future__ import print_function
import msys, sys, math, os

def count_overrides(nb):
    ''' return the number of pairwise nonbonded interactions which 
    are affected by the overrides in the given term table. '''
    # count the number of atoms for each nonbonded param
    nparams = nb.params.nparams
    pcount = [0] * nparams
    for t in nb.terms:
        pcount[t.param.id] += 1
    cnt = sum(pcount[pi.id] * pcount[pj.id] for pi, pj in nb.overrides())
    return cnt
        
def duplicate_overrides(nb, pdict):
    # duplicate the overridden parameters and reassign while disambiguating
    nb.params.addProp('override', int)
    override_base = max(p['override'] for p in nb.params.params)+1
    for p, atoms in list(pdict.items()):
        p_ = p.duplicate()
        p['override']=p.id+override_base
        p_['override']=p_.id+override_base
        for a in atoms: nb.term(a.id).param = p_
        # for the parameters that we duplicated, if there was already an
        # override involving that type, duplicate it with the new type.
        for (pi, pj), op in list(nb.overrides().items()):
            if pi==p: pi=p_
            if pj==p: pj=p_
            nb.setOverride(pi,pj,op)

def add_override(nb, **params):
    p = nb.override_params.addParam()
    for k,v in list(params.items()): 
        nb.override_params.addProp(k,type(v))
        p[k]=v
    return p

def apply_override(nb, sel1, sel2, p):
    t1 = set(nb.term(a.id).param for a in sel1)
    t2 = set(nb.term(a.id).param for a in sel2)
    for pi in t1:
        for pj in t2:
            nb.setOverride(pi,pj,p)

def Override(mol, sel1, sel2, **params):
    ''' override the nonbonded interaction for selections sel1 and sel2
    with the given keyword parameters.  '''

    # add properties
    nb = mol.table('nonbonded')
    pdict = dict()
    for a in sel1+sel2:
        p = nb.term(a.id).param
        pdict.setdefault(p,[]).append(a)

    duplicate_overrides(nb, pdict)
    p = add_override(nb, **params)
    apply_override(nb, sel1, sel2, p)


def combine_arithmetic_geometric(pi,pj):
    sigma =        0.5*(pi['sigma'] + pj['sigma'])
    epsilon = math.sqrt(pi['epsilon'] * pj['epsilon'])
    return sigma, epsilon

def combine_geometric(pi,pj):
    sigma = math.sqrt(pi['sigma'] * pj['sigma'])
    epsilon = math.sqrt(pi['epsilon'] * pj['epsilon'])
    return sigma, epsilon

def Scale(mol, sel1, sel2, scale_sigma, scale_epsilon):
    ''' scale the sigma and epsilon interaction for selections s1 and s2
    with the given keyword parameters.  '''

    rule = mol.nonbonded_info.vdw_rule.lower()
    if rule=='arithmetic/geometric': 
        combine = combine_arithmetic_geometric
    elif rule=='geometric': 
        combine = combine_geometric
    else:
        raise ValueError("Don't understand vdw_rule '%s'" % rule)

    # partition selections by type, and call Override for each pair
    nb = mol.table('nonbonded')

    pdict1 = dict()
    pdict2 = dict()
    for a in sel1:
        p = nb.term(a.id).param
        pdict1.setdefault(p,[]).append(a)
    for a in sel2:
        p = nb.term(a.id).param
        pdict2.setdefault(p,[]).append(a)

    duplicate_overrides(nb, pdict1)
    duplicate_overrides(nb, pdict2)

    for p1, s1 in list(pdict1.items()):
        for p2, s2 in list(pdict2.items()):
            sigma, epsilon = combine(p1,p2)
            p = add_override(nb, 
                             sigma=scale_sigma*sigma, 
                             epsilon=scale_epsilon*epsilon)
            apply_override(nb, s1, s2, p)


scale_doc = '''
dms-scale-vdw input.dms output.dms [ options ]

*dms-scale-vdw* scales the vdw interactions between multiple groups of atoms.
The vdw interactions between each ligand group will be scaled by the
specified amount.  As many ligands may be specified as desired, though
different implementations on Desmond and Anton may in practice limit the
number possible

Currently, the vdw functional form of the DMS file must be "vdw_12_6".  

This tool uses the `nonbonded_combined_param` table in the DMS file to store
the overrides and therefore should not be used with versions of Anton
software older than 2.9.2  
'''

def scale_main():
    import optparse
    parser = optparse.OptionParser(scale_doc)

    parser.add_option('-s', '--scale-sigma', 
            help='scale factor for sigma')
    parser.add_option('-e', '--scale-epsilon',
            help='scale factor for epsilon')
    parser.add_option('-l', '--ligand', action='append', default=[],
            help='atom selection for a ligand')

    opts, args = parser.parse_args()
    if len(args)!=2:
        parser.error("incorrect number of arguments")

    n=len(opts.ligand)
    if n<2:
        parser.error("need at least 2 ligand selections")

    if opts.scale_sigma is None: parser.error("--scale-sigma is required")
    if opts.scale_epsilon is None: parser.error("--scale-epsilon is required")
    sig = float(opts.scale_sigma)
    eps = float(opts.scale_epsilon)

    print("Sigma scale:   %f" % sig)
    print("Epsilon scale: %f" % eps)

    ifile, ofile = args
    print("Loading DMS file <%s>" % ifile)
    mol=msys.LoadDMS(ifile)

    for i in range(0,n):
        s1 = mol.select(opts.ligand[i])
        name = '1st' if i==0 else '2nd' if i==1 else '3rd' if i==2 else '%sth' % (i+1)

        if not s1:
            print("ERROR: No atoms in the %s selection: %s" % (
                    name, opts.ligand[i]))
            exit(1)
        else:
            print("Selected %d atoms with %s selection" % (len(s1), name))
        for j in range(i+1,n):
            s2 = mol.select(opts.ligand[j])
            vdw.Scale(mol, s1, s2, sig, eps)

    mol.coalesceTables()
    mol = mol.clone()
    print("Writing DMS file <%s>" % ofile)
    msys.SaveDMS(mol, ofile)
    print("OK")

override_doc = '''
dms-override-vdw input.dms output.dms [ options ]

Override vdw interactions between selected atoms.

*dms-override-vdw* changes the vdw interaction between two specified groups
of atoms to the specified values of sigma and epsilon.  All options (sigma,
epsilon, selection0, selection1) are required, and the selection groups must
not be empty.  

Currently, the vdw functional form of the DMS file must be "vdw_12_6".  

This tool uses the `nonbonded_combined_param` table in the DMS file to store
the overrides and therefore should not be used with versions of Anton
software older than 2.9.2  
'''

def override_main():
    import optparse
    parser = optparse.OptionParser(__doc__)

    parser.add_option('--sigma', help='Vdw sigma')
    parser.add_option('--epsilon', help='Vdw epsilon')
    parser.add_option('--selection0', help='atom selection for first group')
    parser.add_option('--selection1', help='atom selection for second group')

    opts, args = parser.parse_args()
    if len(args)!=2:
        parser.error("incorrect number of arguments")
    try:
        sigma = float(opts.sigma)
    except TypeError:
        print("Missing/invalid sigma '%s'" % (opts.sigma), file=sys.stderr)
        exit(1)
    try:
        epsilon = float(opts.epsilon)
    except TypeError:
        print("Missing/invalid epsilon '%s'" % (opts.epsilon), file=sys.stderr)
        exit(1)

    ifile, ofile = args
    print("Loading DMS file <%s>" % ifile)
    mol=msys.LoadDMS(ifile)
    s1=str(opts.selection0)
    s2=str(opts.selection1)
    s1=mol.select(s1)
    s2=mol.select(s2)
    if not s1:
        print("ERROR: No atoms in selection0:", opts.selection0)
    else:
        print("Selected %d atoms in selection0" % len(s1))
    if not s2:
        print("ERROR: No atoms in selection1:", opts.selection1)
    else:
        print("Selected %d atoms in selection1" % len(s2))
    if s1 and s2:
        vdw.Override(mol, s1, s2, **{'sigma' : sigma, 'epsilon' : epsilon })
        mol.coalesceTables()
        mol = mol.clone()
    else:
        print("Exiting abnormally")
        exit(1)
    print("Writing DMS file <%s>" % ofile)
    msys.SaveDMS(mol, ofile)
    print("OK")

