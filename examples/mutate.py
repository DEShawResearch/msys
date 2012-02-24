#!/usr/bin/env desres-exec
#{
# desres-cleanenv \
# -m Python/2.7.1-06A/bin \
# -m msys/1.0.10/lib-python \
# -- python $0 "$@"
#}

import sys, msys
from msys import builder

# load a charmm topology file
defs = builder.Defs()
defs.load('top_all27_prot_lipid_na.inp')

# load the input structure
mol=msys.LoadDMS(sys.argv[1], structure_only=True)

# split into protein and everything else
pro=mol.clone('protein')
other=mol.clone('not protein')

# rename residues to match topology definitions
for old, new in {
        'HID'  : 'HSD',
        'LYSH' : 'LYS',
        'CYSH' : 'CYS',
        }.items():
    for a in pro.select('resname %s' % old):
        a.residue.name=new

# rename protein H to HN
for a in pro.select('name H'): a.name='HN'

# mutate residue 7 to an ALA
res=pro.residue(7)
print "mutating %s %d to ALA" % (res.name, res.resid)
res.name='ALA'

# rebuild
new=defs.build(pro)

# append the rest
new.append(other)

# write back out.
msys.SaveDMS(new, 'out.dms')

