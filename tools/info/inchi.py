import msys
import subprocess as SP

def Strings(mol, add_hydrogen=False, include_stereo=False,
                     without_fixed_hydrogen=False):
    ''' invokes the inchi-1 program and returns the inchi strings for
    each unique fragment in the chemical system. '''

    args='inchi-1 -STDIO -AuxNone -NoLabels'.split()

    if not add_hydrogen:
        args.append('-DoNotAddH')
    if not include_stereo:
        args.append('-SNon')
    if not without_fixed_hydrogen:
        args.append('-FixedH')

    p = SP.Popen(args, stdin=SP.PIPE, stdout=SP.PIPE, stderr=SP.PIPE,
            bufsize=4096)
    for fragid in msys.FindDistinctFragments(mol):
        frag = mol.clone('fragid %d' % fragid)
        if frag.natoms > 999:
            continue
        msys.SaveSDF(frag, p.stdin)
    p.stdin.close()
    if p.wait():
        raise RuntimeError, p.stderr.read()
    return p.stdout.read().strip().split('\n')

