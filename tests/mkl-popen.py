#!/usr/bin/garden-exec
#{
### garden env-keep-only
### echo $$
### garden load desres-python/2.7.15-04c7/bin
### exec python $0 "$@"
#}

'''
test for MKL-related crash
'''

from __future__ import print_function
import os

def popen():
    import os
    cmd='ps -p %d h -o vsz' % os.getpid()
    with os.popen(cmd) as p:
        s = p.read()
    return int(s)

def main():
    print("try 1" )
    popen()

    import numpy
    numpy.linalg.svd(numpy.ones((300,300)))

    print("try 2")
    popen()

if __name__=="__main__":
    main()

