import msys

# convert residue name to one-letter code
code = {
        "ALA" : 'A',
        "ASX" : 'B',    # Asn/Asp
        "CYM" : 'C',    # Negative Cys
        "CYS" : 'C',
        "CYX" : 'C',    # S-S bonded Cys
        "ASH" : 'D',    # neutral Asp
        "ASP" : 'D',
        "GLH" : 'E',    # neutral Glu
        "GLU" : 'E',
        "PHE" : 'F',
        "GLY" : 'G',
        "HSD" : 'H',    # Delta His
        "HSE" : 'H',    # Epsilon His
        "HID" : 'H',    # Delta His
        "HIE" : 'H',    # Epsilon His
        "HIP" : 'H',    # Positive His
        "HIS" : 'H',
        "ILE" : 'I',
        "LYN" : 'K',    # neutral Lys
        "LYS" : 'K',
        "LEU" : 'L',
        "MET" : 'M',
        "ASN" : 'N',
        "PRO" : 'P',
        "GLN" : 'Q',
        "ARG" : 'R',
        "SER" : 'S',
        "THR" : 'T',
        "VAL" : 'V',
        "TRP" : 'W',
        "UNK" : 'X',
        "TYM" : 'Y',    # neutral Tyr
        "TYR" : 'Y',
        "GLX" : 'Z',    # Gln/Glu
        }

def Sequences(system_or_chain, distinct=True):
    ''' return list of sequences, one for each chain, for the given input.
    The sequence will be returned as a string, with characters corresponding
    to the 1-letter sequence codes.  If a Chain is provided instead of a
    System, only one sequence will be returned.

    If distinct is True, only distinct sequences will be returned.
    '''
    if isinstance(system_or_chain, msys.Chain):
        chains = [system_or_chain]
    elif isinstance(system_or_chain, msys.System):
        chains = system_or_chain.chains
    else:
        t=type(system_or_chain)
        raise TypeError, "Expected System or Chain, got %s" % t
    seqs = []
    for c in chains:
        seq = ''.join(code.get(r.name,'X') for r in c.residues)
        seq = seq.strip('X')
        if not seq.translate(None, 'X'):
            # got all X
            seq = ''
        seqs.append(seq)
    if distinct:
        seqs = sorted(set(seqs))
    return seqs


