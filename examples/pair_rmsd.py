#!/usr/bin/garden-exec
#{
# exec garden with -c -m msys/1.7.117/lib-python -- python $0 "$@"
#}


import msys, sys

def main():
    f1, f2 = sys.argv[1:]
    for mol1, mol2 in zip(msys.LoadMany(f1), msys.LoadMany(f2)):
        print mol1.atomsel('noh').alignedRMSD(mol2.atomsel('noh'))

main()



