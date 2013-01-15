#!/usr/bin/env desres-exec
#{
# desres-cleanenv -m periodicfix/2.4.7/lib-python -- python $0 "$@"
#}

import os, sys, unittest
TMPDIR=os.getenv('TMPDIR', 'objs/Linux/x86_64')
sys.path.insert(0,os.path.join(TMPDIR, 'lib', 'python'))

import msys, sys
import periodicfix as PF
import numpy as NP

def compute_aligned_rmsd(rpos, ratoms, tpos, patoms):
    rx = rpos[[a.id for a in ratoms]]
    px = tpos[[a.id for a in patoms]]
    rx -= rx.mean(0)
    px -= px.mean(0)
    mat = PF.compute_alignment(rx, px, None)
    tx = NP.dot(px, mat)
    tx -= rx    # deltas
    tx *= tx    # square deltas
    return tx.sum()

def fix_atom_order(ref, tgt, refsel, targetsel):

    rpos = ref.positions.astype('f')
    tpos = tgt.positions.astype('f')
    # find the set of residues in each selection
    refres = set(a.residue.id for a in ref.select(refsel))
    tgtres = set(a.residue.id for a in tgt.select(targetsel))
    assert len(refres)==len(tgtres), "System selections contain different numbers of residues: ref=%d, target=%d" % (len(refres), len(tgtres))

    # create a mapping from target residues to reference residues.  We are
    # assuming that the order of residue ids is the same in both cases.
    mapping = dict(zip(sorted(tgtres), sorted(refres)))

    # iterate over target residues.  If there is no analog in refres,
    # use the existing ordering; otherwise do a graph isomorphism.
    atoms=[]
    for t_res in tgt.residues:
        tatoms = t_res.atoms
        r = mapping.get(t_res.id)
        if r is not None:
            r_res = ref.residue(r)
            ratoms = r_res.atoms
            assert len(ratoms)==len(tatoms), "Residue %s in target has different number of atoms (%d) than residue %s in reference (%d)" % (t_res, len(tatoms), r_res, len(ratoms))
            rg = msys.Graph(ratoms)
            tg = msys.Graph(tatoms)

            perms = []
            for d in rg.matchAll(tg):
                # permuted atoms in d
                patoms = [d[r] for r in ratoms]
                # aligned distance
                dist = compute_aligned_rmsd(rpos, ratoms, tpos, patoms)
                perms.append((dist, patoms))

            if not perms:
                print "ref atoms:"
                for a in sorted(ratoms, key=lambda x: x.name):
                    print a.name, a.atomic_number, a.nbonds
                print "tgt atoms:"
                for a in sorted(tatoms, key=lambda x: x.name):
                    print a.name, a.atomic_number, a.nbonds
                raise RuntimeError, "No isomorphism for residue %s in ref and residue %s in target" % (
                        r_res, t_res)

            perms.sort()
            atoms.extend(perms[0][1])
        else:
            atoms.extend(t_res.atoms)
    return tgt.permuted(atoms)

def main():
    import optparse
    parser = optparse.OptionParser(__doc__)
    parser.add_option('--refsel', default='all', 
            help='selected atoms in reference structure')
    parser.add_option('--targetsel', default='all', 
            help='selected atoms in target structure')

    opts, args = parser.parse_args()
    if len(args)!=3:
        parser.error("incorrect number of arguments")

    rfile, ifile, ofile = args
    
    ref = msys.Load(rfile)
    tgt = msys.Load(ifile)
    out = fix_atom_order(ref, tgt, **opts.__dict__)
    msys.SaveDMS(out, ofile)

if __name__=="__main__": main()


