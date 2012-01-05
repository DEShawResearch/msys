#!/usr/bin/env desres-exec
#{
# desres-cleanenv \
# -m Python/2.7.1-06A/bin \
# -m ndiff/2.00-02/bin \
# -m sqlite/3.7.8.0-03A/bin \
# -- python $0 "$@"
#}

'''
read dms files in, write them back out, then compare using
dms-diff.
'''

import sys, os
sys.path.insert(0,'objs/Linux/x86_64/lib/python')
import msys

DUMP='./dmsdump.py --without-provenance --without-forcefield --without-paraminfo'
DIFF='ndiff -quiet -relerr 1e-4'

import multiprocessing as MP

def doit(path):
    base=os.path.basename(path)
    try:
        print base
        input=path+'-in'
        output=path+'-out'
        os.system('cp %s %s' % (path, input))
        os.system('chmod u+w %s' % input)
        os.system('sqlite3 %s "update particle set name=trim(name)"' % input)

        mol=msys.LoadDMS(input)
        msys.SaveDMS(mol,output)


        os.system('%s %s > %s.dump' % (DUMP, input, input))
        os.system('%s %s > %s.dump' % (DUMP, output, output))

        rc=os.system('%s %s.dump %s.dump' % (DIFF, input, output))
        return base, rc 
    except RuntimeError, e:
        return base, e

if __name__=="__main__":
    p=MP.Pool(8)
    results = p.map(doit, sys.argv[1:])
    failed = [r for r in results if r[1] is not 0]
    print "---------------------------------------------------------------"
    print "%d/%d files passed." % (len(results)-len(failed), len(results))
    if failed:
        for path, rc in failed:
            print "%s: %s" % (path, rc)
        exit(1)

