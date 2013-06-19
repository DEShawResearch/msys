import _atommatch
import sys
import msys
import numpy as np

def default_score_fct(atoms1, atoms2):
    ''' Default scoring function for use in AtomMatch.

    Rule 1: Do not match hydrogens to other elements.
    Rule 2: Do not match if any bond lengths differ by more than 10%.
    Rule 3: Otherwise, score is the number of matched atoms.
    '''
    from collections import deque
    # Rule 1: Do not match hydrogens to other elements
    for a, b in zip(atoms1, atoms2):
        if a.atomic_number == 1 and b.atomic_number != 1 \
                or a.atomic_number != 1 and b.atomic_number == 1:
            return -np.inf
    # Rule 2: Do not match if bond lengths differ by more than 10%
    for i in range(len(atoms1)-1):
        bonded = atoms1[i].bonded_atoms
        for j in range(i+1, len(atoms1)):
            if atoms1[j] in bonded:
                dist1 = np.sqrt(sum((atoms1[i].pos - atoms1[j].pos)**2))
                dist2 = np.sqrt(sum((atoms2[i].pos - atoms2[j].pos)**2))
                if abs(dist1 - dist2) > 0.1 * min(dist1, dist2):
                    return -np.inf
    # Rule 3: Otherwise score is number of matched atoms
    return len(atoms1)

# Computes RMSD alignment
def _rawalign(X, Y, w=None):
    from numpy.linalg.linalg import svd, det
    if w == None:
        w = np.ones(X.shape[0], np.float32)
    w *= (1.0 / w.sum())
    xm = np.inner(X.T, w)
    ym = np.inner(Y.T, w)
    Xo = X - xm # N x 3, centered x
    Yo = Y - ym # N x 3, centered y
    u, s, v = svd(np.inner(Yo.T, Xo.T * w)) # Law: dot(u, dot(diag(s), v)) = a
    trxxtw = np.inner((Xo * Xo).sum(axis=1), w) # tr(XX^TW)
    tryytw = np.inner((Yo * Yo).sum(axis=1), w) # tr(YY^TW)
    Ap = np.array([1, 1, det(np.dot(u, v))]) # diagonal values
    sqerr = trxxtw + tryytw - 2 * np.inner(Ap, s) # tr(A_p \Sigma)
    rmsd = np.sqrt(max(sqerr, 0))
    A = np.dot(v.T, (u * Ap).T)
    b = ym - np.dot(xm, A)
    return A, b, rmsd

def AtomMatch(mol1, mol2, sel1='all', sel2='all', score_fct=default_score_fct):
    '''Get best mapping of atoms between two molecules.

    Returns the best isomorphism between subgraphs in mol1 and mol2,
    where "best" is determined by the following rules:

    (1) Each subgraph is a union of complete biconnected components of
        the molecule (so we cannot match only part of a ring system).

    (2) The match has the highest score among all possible subgraph
        matches satisfying (1), where the score of a match is the sum
        of scores over matching biconnected components as determined by
        'score_fct'. The null match has score 0.

    (3) Among multiple matches satisfying (2), the match contains the
        most matched atoms.

    (4) Among multiple matches satisfying (3), the match minimizes RMSD
        between the matched atoms.

    Rules (1)-(3) are enforced exactly. As the number of matches that
    satisfy rules (1)-(3) may increase exponentially in the number of
    symmetries of the molecule, rule (4) is not enforced exactly but
    rather is implemented by a heuristic search.

    Rule (2) is determined by a user-defined score function of the form

        score_fct(atoms1, atoms2) -> float

    where atoms1 and atoms2 are isomorphic lists of msys.Atom's from
    mol1 and mol2 from a single biconnected component of each molecule,
    and the return value is a numeric score. Score functions may use
    any atom or bond properties available in msys.

    Arguments:
        mol1 -- msys.System
        mol2 -- msys.System
        sel1 -- str, atom selection for mol1
        sel2 -- str, atom selection for mol2
        score_fct --
          f([Atom, ..., Atom], [Atom, ..., Atom]) -> float

    Returns: (match, aligned_mol)

        match -- [(Atom, Atom), ... ], where the first atom in each
        pair is an atom from mol1 and the second is the matching atom
        from aligned_mol

        aligned_mol -- msys.System, a clone of the selection sel2 of mol2
        that is RMSD-aligned to mol1 using the positions of the matched
        atoms

    '''
    # Clone selections
    atoms1 = mol1.select(sel1)
    atoms2 = mol2.select(sel2)
    clone1 = msys.CloneSystem(atoms1)
    clone2 = msys.CloneSystem(atoms2)
    if len(clone1.updateFragids()) != 1 or len(clone2.updateFragids()) != 1 \
            or len(set.union(*[set(a.bonded_atoms)
                for a in atoms1])) != len(atoms1) \
            or len(set.union(*[set(a.bonded_atoms)
                for a in atoms2])) != len(atoms2):
        raise RuntimeError, \
                'Atom selection must correspond to a single complete fragment'
    # Call C++ function
    def f(m1, a1, m2, a2):
        atoms1 = [msys.Atom(m1, id) for id in a1]
        atoms2 = [msys.Atom(m2, id) for id in a2]
        return score_fct(atoms1, atoms2)
    rep = _atommatch.ScoreFct(f)
    out = _atommatch.AtomMatch(clone1._ptr, clone2._ptr, rep)
    match = [(atoms1[id1], clone2.atom(id2)) for (id1, id2) in out]
    # Align mol2 to mol1 based on matched atoms
    mol2_copy = mol2.clone()
    A, b, rmsd = _rawalign(np.array([a2.pos for a1,a2 in match]),
            np.array([a1.pos for a1,a2 in match]))
    for a2 in clone2.atoms:
        a2.pos = np.dot(a2.pos, A) + b
    return match, clone2
