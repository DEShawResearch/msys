#!/usr/bin/env desres-exec
#{
# desres-cleanenv \
# -m Python/2.7.1-06A/bin \
# -m msys/1.0.8/lib-python \
# -- python $0 "$@"
#}

import msys, math

def check(mol):
    for b in mol.bonds:
        ai=b.first
        aj=b.second
        ipos=ai.pos
        jpos=aj.pos
        r=math.sqrt(sum((xi-xj)**2 for xi, xj in zip(ipos, jpos)))
        if r<1.1 and ai.atomic_number!=1 and aj.atomic_number!=1:
            print "SUSPICIOUSLY SHORT BOND OF LENGTH %s between atoms %d and %d" % (
                    r, ai.id, aj.id)
        if r>2.5:
            print "SUSPICIOUSLY LONG BOND OF LENGTH %s between atoms %d and %d" % (
                    r, ai.id, aj.id)

        
if __name__=="__main__":
    import sys
    for p in sys.argv[1:]:
        mol=msys.Load(p)
        check(mol)
        




