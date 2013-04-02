import sys
import msys
import subprocess as SP
import collections
import os

if os.getenv('DESRES_ROOT'):
    _inchi_exe = 'desres-cleanenv -m inchi/1.04-04 -- inchi-1'
else:
    _inchi_exe = 'inchi-1'
_inchi_exe = _inchi_exe.split()

def Strings(mol, add_hydrogen=False, include_stereo=False,
                     without_fixed_hydrogen=False):
    ''' invokes the inchi-1 program and returns the inchi strings for
    each unique fragment in the chemical system. '''

    args = _inchi_exe + '-STDIO -AuxNone -NoLabels'.split()

    if not add_hydrogen:
        args.append('-DoNotAddH')
    if not include_stereo:
        args.append('-SNon')
    if not without_fixed_hydrogen:
        args.append('-FixedH')

    p = SP.Popen(args, stdin=SP.PIPE, stdout=SP.PIPE, stderr=SP.PIPE,
            bufsize=4096)
    for fragid in msys.FindDistinctFragments(mol):
        frag = mol.clone('fragid %d and atomicnumber>0' % fragid)
        if frag.natoms > 999:
            continue
        msys.SaveSDF(frag, p.stdin)
    p.stdin.close()
    if p.wait():
        raise RuntimeError, p.stderr.read()
    return p.stdout.read().strip().split('\n')

StandardInchi = collections.namedtuple(
    'StandardInchi', 'InChI AuxInfo InChIKey'
    )

def InChI(mol, FixedH=False, SNon=False, DoNotAddH=False):
    ''' compute the inchi for the given component and return
    the InChI, InChIKey, and AuxInfo as a named tuple.  
    '''
    args = _inchi_exe + '-KEY -STDIO -NoLabels'.split()
    if FixedH: args.append('-FixedH')
    if SNon: args.append('-SNon')
    if DoNotAddH: args.append('DoNotAddH')

    p = SP.Popen(args, stdin=SP.PIPE, stdout=SP.PIPE, stderr=SP.PIPE,
                        bufsize=4096)
    msys.SaveSDF(mol, p.stdin)
    p.stdin.close()
    rc = p.wait()
    log = p.stderr.read()
    if rc:
        raise RuntimeError, log
    output = p.stdout.read()
    lines = output.strip().split('\n')
    if len(lines)!=3:
        raise RuntimeError, "Got %d lines instead of 3: %s" % (len(lines), output)
    inchi, aux, key = lines
    return StandardInchi(
            InChI=inchi[6:], 
            AuxInfo=aux[8:],
            InChIKey=key[9:],
            )

