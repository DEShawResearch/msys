
import msys, sys
from collections import Counter

SP = msys.SmartsPattern


amino='[$([NX3H2,NX4H3+]),$([NX3H](C)(C))][CX4H](*)[CX3](=[OX1])[OX2H,OX1-,N]'

sides={
        'A' : '[CH3X4]',
        'C' : '[CH2X4][SX2H,SX1H0-]',
        'D' : '[CH2X4][CX3](=[OX1])[OH0-,OH]',
        'E' : '[CH2X4][CH2X4][CX3](=[OX1])[OH0-,OH]',
        'F' : '[CH2X4][cX3]1[cX3H][cX3H][cX3H][cX3H][cX3H]1',
        'H' : '[CH2X4][#6X3]1:[$([#7X3H+,#7X2H0+0]:[#6X3H]:[#7X3H]),$([#7X3H])]:[#6X3H]:[$([#7X3H+,#7X2H0+0]:[#6X3H]:[#7X3H]),$([#7X3H])]:[#6X3H]1',
        'I' : '[CHX4]([CH3X4])[CH2X4][CH3X4]',
        'K' : '[CH2X4][CH2X4][CH2X4][CH2X4][NX4+,NX3+0]',
        'L' : '[CH2X4][CHX4]([CH3X4])[CH3X4]',
        'M' : '[CH2X4][CH2X4][SX2][CH3X4]',
        'N' : '[CH2X4][CX3](=[OX1])[NX3H2]',
        #'Q' : '[CH2X4][CH2X4][CX3](=[OX1])[OH0-,OH]',
        'R' : '[CH2X4][CH2X4][CH2X4][NHX3][CH0X3](=[NH2X3+,NHX2+0])[NH2X3]',
        'S' : '[CH2X4][OX2H]',
        'T' : '[CHX4]([CH3X4])[OX2H]',
        'V' : '[CHX4]([CH3X4])[CH3X4]',
        'W' : '[CH2X4][cX3]1[cX3H][nX3H][cX3]2[cX3H][cX3H][cX3H][cX3H][cX3]12',
        'Y' : '[CH2X4][cX3]1[cX3H][cX3H][cX3]([OHX2,OH0X1-])[cX3H][cX3H]1',
        }

gly='[$([$([NX3H2,NX4H3+]),$([NX3H](C)(C))][CX4H2][CX3](=[OX1])[OX2H,OX1-,N])]'
pro='[$([NX3H,NX4H2+]),$([NX3](C)(C)(C))]1[CX4H]([CH2][CH2][CH2]1)[CX3](=[OX1])[OX2H,OX1-,N]'

aaname={
        'A' : 'ALA',
        'C' : 'CYS',
        'D' : 'ASP',
        'E' : 'GLU',
        'F' : 'PHE',
        'G' : 'GLY',
        'H' : 'HIS',
        'I' : 'ILE',
        'K' : 'LYS',
        'L' : 'LEU',
        'M' : 'MET',
        'N' : 'ASN',
        'P' : 'PRO',
        'Q' : 'GLN',
        'R' : 'ARG',
        'S' : 'SER',
        'T' : 'THR',
        'V' : 'VAL',
        'W' : 'TRP',
        'Y' : 'TYR',
        }




# generate smarts patterns
smarts=dict()
smarts['G']=SP(gly)
smarts['P']=SP(pro)
for k,v in sides.iteritems():
    smarts[k] = SP(amino.replace('*', v))


for p in sys.argv[1:]:
    mol=msys.Load(p).clone('not water and not ion and not lipid')
    print "Loaded", mol.name
    SP.Annotate(mol)

    # make a lookup from atom to residue
    rmap = [a.residue.id for a in mol.atoms]

    # apply each pattern, generating a mapping from atomid to code.
    atommap=dict()
    for k,p in smarts.iteritems():
        for m in p.findMatches(mol):
            # sometimes these matches are spanning multiple residues.
            # hack around this by partitioning the match by residue and
            # keeping only the largest group.
            c = Counter(rmap[i] for i in m).most_common()
            c = [e for e in c if e[1] == c[0][1]]
            if len(c)==1:
                r = c[0][0]
            else:
                # no unique largest group
                print "ERROR:", c
            for i in m:
                if rmap[i]==r:
                    atommap.setdefault(i,[]).append(k)

    for r in mol.residues:
        # Find the pattern with the longest hit on each residue
        ids = [a.id for a in r.atoms]
        m = Counter(sum((atommap.get(i,[]) for i in ids), [])).most_common()
        if not m:
            code = 'X'
        else:
            m = [e for e in m if e[1]==m[0][1]]
            if len(m)==1:
                code = m[0][0]
            else:
                print "ERROR: residue %s %d" % (r.name, r.resid), m

        if r.name[:2] != aaname.get(code, 'X')[:2]:
            print "MISMATCH: %5d %3s != %s" % ( r.resid,r.name,aaname.get(code))
    
