#!/usr/bin/env desres-exec
#{
# desres-cleanenv \
# -m Python/2.7.1-06A/bin \
# -m msys/1.0.10/lib-python \
# -- python $0 "$@"
#}

import msys
from msys import builder, solvate, neutralize

mol=msys.Load('dodecyl_beta_maltoside.pdb')

# put all the atoms into the same residue
for a in mol.atoms: a.residue = mol.residue(0)

# Load the topology file
defs = builder.Defs()
defs.load('maltoside-waterC36.top')

# Build everything in the system.
built=defs.build(mol)

# replicate 5x5x5
replicas=msys.CreateSystem()
delta = 20.0    # spacing between replicas
resid=1
for i in range(5):
    for j in range(5):
        for k in range(5):
            rep=built.clone()
            rep.residue(0).resid = resid
            resid += 1
            rep.translate((delta*i, delta*j, delta*k))
            replicas.append(rep)

# recenter on the origin
replicas.translate(-x for x in replicas.center)

# solvate; this also sets global cell size
solvated = solvate.Solvate(replicas, dims=[5*delta], verbose=True)

# neutralize
neutralize.Neutralize(solvated, concentration=0.050, verbose=True)

# write back out.
msys.SaveDMS(solvated, 'out.dms')

