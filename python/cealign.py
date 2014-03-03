
from msys._cealign import CEAlign
import itertools

def gen_matches(A,B):
    ''' generate all possible ways of matching elements of A with elements of
    B, using as many members of the smaller set as possible.  For example,
    if A={a,b} and B={x,y,z}, the possible matches are:
        ax-by ax-bz ay-bx ay-bz az-bx az-by
    '''
    n=max(len(A), len(B))
    k=min(len(A), len(B))
    for p in itertools.permutations(range(n), k):
        if len(A)>len(B):
            # take p elements from A and zip with B
            yield zip((A[i] for i in p), B)
        else:
            # take p elements from B and zip with A
            yield zip(A, (B[i] for i in p))


def cealign(ref, target, refsel='name CA', targetsel='name CA'):
    ''' Use the CEAlign algorithm to find a good alignment of the atoms
    in target onto the atoms in ref.  

    If ref or target comprise multiple disconnected chains, all possible
    matches are considered.
    '''
    rmol = ref.clone(refsel)
    tmol = target.clone(targetsel)
    rpos = rmol.getPositions()
    tpos = tmol.getPositions()

    # chainhash stores the result of calling cealign on a single
    # pair of chains.  The values are lists of paths, and each path
    # is a pair of atom ids lists.  I.e. (A1,A2) -> [(r1,t1),(r2,t2), ...]
    chainhash = dict()
    engine = CEAlign.WithDefaults()
    pathlists=[]

    # for each matching of chains in ref to chains in target...
    for chain_pairs in gen_matches(rmol.chains, tmol.chains):
        mapping=[]
        for chain_pair in chain_pairs:
            pathlist = chainhash.get(chain_pair)
            if pathlist is None:
                rchain, tchain = chain_pair
                ratoms = [a.id for a in (r for r in rchain.residues)]
                tatoms = [a.id for a in (r for r in tchain.residues)]
                pathlist = engine.compute(ratoms, rpos, tatoms, tpos)
                chainhash[chain_pair] = pathlist
            mapping.append(pathlist)
        # by construction, all paths in a given pathlist are of the same
        # length, so the elements of the cartesian product have the same
        # length, too.  
        print mapping

