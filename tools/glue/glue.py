
import msys

def _find_glue(sys, ids):
    ''' Return minimal glue for the given ids '''
    frags_outside = set(sys.atom(i).fragid for i in ids)
    if not frags_outside: return []
    frags_inside = set()
    latest_frag = frags_outside.pop()
    pairs = []
    sel='index ' + ' '.join(str(i) for i in ids)

    while frags_outside:
        frags_inside.add(latest_frag)
        inside_frag_string = ' '.join('%d' % s for s in frags_inside)
        outside_frag_string = ' '.join('%d' % s for s in frags_outside)
        r=1
        done=False

        while not done:
            query='(%s) and (fragment %s) and (within %s of fragment %s)' % (
                    sel, outside_frag_string, r, inside_frag_string)
            #print "query:", query
            indices=sys.select(query)
            #print "indices:", [x for x in indices]
            if not indices:
                r += 1
                continue
            # Found the next fragment to add to the MST.
            # Find which atom in frags_inside is closest to outside_index
            outside_index = indices[0]
            latest_frag = sys.atom(outside_index).fragid
            ox = sys.atom(outside_index).x
            oy = sys.atom(outside_index).y
            oz = sys.atom(outside_index).z
            assert latest_frag not in frags_inside
            frags_outside.remove(latest_frag)
            best_d = None
            best_atom = None
            for inside_index in sys.select('(%s) and fragment %s' % (
                sel, inside_frag_string)):
                ix = sys.atom(inside_index).x
                iy = sys.atom(inside_index).y
                iz = sys.atom(inside_index).z
                d = (ox-ix)**2 + (oy-iy)**2 + (oz-iz)**2
                #print "outside %d inside %d d %f\n" % (
                        #outside_index, inside_index, d)
                if best_atom is None or d < best_d:
                    best_atom = inside_index
                    best_d = d
            assert best_atom is not None
            pairs.append((best_atom, outside_index))
            done=True
            break

    return pairs

def FindGlue(atoms):
    ''' Return ids of the atoms forming minimal glue for the given input
    set of atoms '''
    sys, ids = msys._find_ids(atoms)
    return _find_glue(sys, ids)

