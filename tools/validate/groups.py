
import numpy

def check_midpoint(table, radius):
    print "Checking midpoints for table %s with %d atoms, %d terms" % (
        table.name, table.natoms, table.nterms)
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
                print "  Term %d too large with atoms %s" % (
                t, ' '.join(map(str, atoms)))
                break
    if bad:
        print "Affected atoms:", ' '.join(map(str, sorted(bad)))
    return not bad

def check_replicate(table, radius):
    return True

def check(table, radius):
    if table.category in ('bond', 'exclusion'):
        return check_midpoint(table, radius)
    elif table.category in ('constraint', 'virtual', 'polar'):
        return check_replicate(table, radius)
    return False

