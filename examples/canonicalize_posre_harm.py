#!/usr/bin/env desres-exec
#{
# desres-cleanenv \
# -m Python/2.7.1-06A/bin \
# -m msys/1.3.7/lib-python \
# -- python $0 "$@"
#}

'''
Some posre_harm tables have x0, y0, z0 as param properties rather than
term properties.  This script enforces the convention that these properties
ought to be term properties.
'''

import msys

def canonicalize(mol):
    if 'posre_harm' not in mol.table_names: 
        print "No posre_harm table found"
        return

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
