#!/usr/bin/env python2.7

import sys, os
sys.path.insert(0,'objs/Linux/x86_64/lib/python')
import msys

import multiprocessing as MP

def doit(path):
    base=os.path.basename(path)
    try:
        print base
        msys.LoadMAE(path)
        return base, None
    except RuntimeError, e:
        return base, e
    except:
        return base, "Unknown error"

if __name__=="__main__":
    p=MP.Pool(8)
    results = p.map(doit, sys.argv[1:])
    failed = [r for r in results if r[1] is not None]
    print "---------------------------------------------------------------"
    print "%d/%d files passed." % (len(results)-len(failed), len(results))
    if failed:
        for path, rc in failed:
            print "%s: %s" % (path, rc)
        exit(1)

