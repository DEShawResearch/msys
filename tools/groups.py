
import numpy

def midpoint_violations(table, radius):
    pos = table.system.getPositions()
    center = numpy.zeros(3, 'd')
    scale = 1.0/table.natoms
    r2 = radius * radius
    bad = set()
    for t in table._ptr.terms():
        center[:] = 0
        atoms = table._ptr.atoms(t)
        for a in atoms:
            center += pos[a]
        center *= scale
        for a in atoms:
            delta = pos[a] - center
            rad2 = numpy.dot(delta, delta)
            if rad2 > r2:
                bad.update(atoms)
                break
    return sorted(bad)

def clone_buffer_violations(table, radius):
    ''' return atoms which violate Desmond's clone buffer restriction 
    at the given radius for terms in the given table.
    '''
    if table.natoms < 2: 
        return []
    elif table.category in ('bond', 'exclusion'):
        return midpoint_violations(table, radius)
    return []

