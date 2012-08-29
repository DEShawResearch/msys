#!/usr/bin/env desres-exec
#{
# exec desres-cleanenv \
# -m Python/2.7.1-06A/bin \
# -m numpy/1.5.1-29A/lib-python \
# -m msys/1.6.0/lib-python \
# -- python $0 "$@"
#}

import msys
import math
import string
import sys

'''
Produce a coverage matrix: rows are features, columns are the 
input files.  Add new features by defining a function starting with 
"feature_" at global scope.
'''

def feature_multiple_virtuals_per_host(mol):
    ''' returns True if multiple virtual sites reference the same
    parent atoms.  '''
    parents=set()
    for T in mol.tables:
        if T.category=='virtual':
            for t in T.terms:
                oldsize = len(parents)
                parents.update(a.id for a in t.atoms[1:])
                newsize = len(parents)
                if oldsize != newsize:
                    return True
    return False

def feature_big(mol):
    return mol.natoms > 100000

def custom_features():
    ''' yield name, predicate for custom features '''
    for name, func in globals().iteritems():
        if name.startswith('feature_') and callable(func):
            yield name[8:], func

def gen_features(mol):
    ''' generator of features found in the system '''
    for t in mol.tables:
        yield t.name
    for name, func in custom_features():
        if func(mol):
            yield name

def output(names, results):
    shortnames = dict()
    codes=string.letters
    for i,n in enumerate(names):
        shortnames[n] = codes[i]

    length = max(len(f) for f in results.iterkeys())
    fmt = "%" + str(length+1) + "s"
    for f in sorted(results.iterkeys()):
        r = results[f]
        f = f.lower() if r else f.upper()
        print fmt % f,
        print "".join(shortnames[x] for x in r)

    print ""
    for n in names:
        print shortnames[n], n


def main():
    import sys
    # mapping from feature to list of molecules
    results=dict()
    for name,_ in custom_features(): results[name]=[]
    paths=sys.argv[1:]
    for arg in paths:
        try:
            mol = msys.Load(arg)
        except:
            print >> sys.stderr, "Failed to load %s" % arg
            continue
        for f in gen_features(mol):
            results.setdefault(f,[]).append(arg)
    output(paths,results)

if __name__=="__main__": main()

