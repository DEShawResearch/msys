import os, msys
import numpy

WATRAD = 2.4
WATSEL = 'oxygen'
WATCON = 1.0
WATBOX = os.path.join(os.path.dirname(__file__), '../../../share/data/tip3p.dms')

def remove_periodic_contacts(_mol, npro, dist):
    mol = _mol.clone()
    natoms = mol.natoms
    wat = _mol.clone('index >= %d' % npro)
    mol.append(wat)
    sel='index < %(natoms)d and index >= %(npro)d and within %(dist)f of (index >= %(natoms)d)' % locals()
    pos = mol.positions
    a, b, c = mol.cell
    bad = []
    for i in (-1,0,1):
        for j in (-1,0,1):
            for k in (-1,0,1):
                if i==0 and i==j and i==k: continue
                delta = i*a + j*b + k*c
                pos[natoms:] += delta
                mol.setPositions(pos)
                ids = mol.selectIds(sel)
                bad.extend(ids)
                pos[natoms:] -= delta

    if bad:
        _mol=_mol.clone('not same fragid as index ' + ' '.join(map(str, bad)))
    return _mol

def solvate(solute, solvent, dims,
            min_solvent_dist=WATCON,
            min_solute_dist=WATRAD,
            solvent_selection=WATSEL,
            center_selection='all',
            verbose=False):
    ''' 
    Tile solvent box around solute, removing overlaps.
    Args:
        solute (msys.System): system to be solvated
        solvent (msys.System): system containing solvent molecules
        dims (list[float]): dimensions of desired solvated system
        min_solvent_dist (float): minimum distance between periodic contacts
                                  of solvent atoms
        min_solute_dist (float): minimium distance between solvent and solute
        solvent_selection (str): solvent selection used for checking contacts

    Returns:
        copy of solute containing copies of the solute, and the global cell set to dims.

    The solvent system must have its box size set appropriately so that it can be tiled
    without introducing overlapping atoms.
    '''
    mol = solute.clone()
    wat = solvent.clone()
    npro = mol.natoms

    if center_selection != 'none': 
        ids = mol.selectIds(center_selection)
        if not ids:
            raise ValueError("Center selection '%s' selected no atoms" % center_selection)
        pos = mol.positions
        pos -= pos[ids].mean(0)
        mol.positions = pos

    # put all the water in one ct
    ct=mol.addCt()
    ct.name='solvate'

    dims = numpy.array([float(x) for x in dims])
    if len(dims) != 3:
        raise ValueError("Dims must be given as list of three values")
    mol.cell[:] = numpy.diag(dims)

    watsize = numpy.diag(wat.cell)
    if (watsize <= 0).any():
        raise ValueError("Invalid solvent global cell: %s" % wat.cell)

    center = mol.center
    dmin = center - 0.5*dims
    dmax = center + 0.5*dims
    nrep = (dims / watsize).astype('i') + 1
    shift = -0.5 * (nrep-1)*watsize
    nx, ny, nz = nrep
    watpos = wat.positions
    newpos = [mol.positions]
    
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                delta = shift + watsize * (i,j,k)
                newpos.append(watpos + delta)
                mol._ptr.append(wat._ptr, ct.id)

    mol.positions = numpy.concatenate(newpos)
    mol.updateFragids()

    toonear='pbwithin %s of index < %s' % (min_solute_dist, npro)
    mol=mol.clone(
        'not same fragid as (index >= %d and (%s) and (%s))' % (
            npro, solvent_selection, toonear))
    if verbose: print("After removing overlap, %d solvent atoms" % (
            mol.natoms - npro))

    # remove molecules whose center is outside the desired box
    xmin, ymin, zmin = dmin
    xmax, ymax, zmax = dmax
    hits = mol.select('index >= %d and (x<%s or y<%s or z<%s or x>%s or y>%s or z>%s)' % (npro, xmin,ymin,zmin, xmax,ymax,zmax))
    frags = set(a.fragid for a in hits)
    if frags:
        fragmap = dict()
        sel = '(%s) and fragid %s' % (solvent_selection, ' '.join(map(str,frags)))
        for a in mol.select(sel): fragmap.setdefault(a.fragid, []).append(a.id)
        outside = []
        pos = mol.getPositions()
        pmin = numpy.array((xmin,ymin,zmin))
        pmax = numpy.array((xmax,ymax,zmax))
        for fragid, ids in fragmap.items():
            c = pos[ids].mean(0)
            if (c<pmin).any() or (c>pmax).any(): outside.append(fragid)
        if outside:
            mol = mol.clone('not fragid ' + ' '.join(map(str,outside)))
            if verbose: print("After removing outside solvent molecules, %d solvent atoms" % (mol.natoms - npro))

    # remove overlap with periodic images
    mol = remove_periodic_contacts(mol, npro, min_solvent_dist)
    if verbose: print("after removing periodic clashes, %d solvent atoms" % (mol.natoms-npro))

    # assign the water chain name and water resids
    for i, c in enumerate(mol.ct(ct.id).chains):
        c.name = 'W%d' % (i+1)
        for j,r in enumerate(c.residues):
            r.resid = j+1

    return mol

def parse_dims(args, mol):
    if args.box is not None:
        dims = [float(x) for x in args.box.split(',')]
    else:
        pos = mol.positions
        extent = (pos.max(0)-pos.min(0)).max()
        dims = [extent + 2*args.thickness]*3
    return dims


def parser():
    import argparse
    import os
    default_solvent = WATBOX
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('solute',
            help="Path to chemical system to be solvated")
    parser.add_argument('-s', '--solvent', default=os.path.realpath(default_solvent),
            help="Path to chemical system containing solvent molecules")
    parser.add_argument('-o', '--output', default='solvate.dms',
            help="Path to output file")
    parser.add_argument('-c', '--center-selection', default='all',
            help="Center solute using given selection; 'none' to skip")
    parser.add_argument('--structure-only', action='store_true',
            help="Load solute and solvent without any forcefield information")

    grp = parser.add_mutually_exclusive_group()
    grp.add_argument('-b', '--box',
            help='water box dimensions: 1 or 3 comma-separated values')
    grp.add_argument('-t', '--thickness', type=float, default=5.0,
            help='Minimum distance between solute and edge of water box')

    parser.add_argument('--solute-dist', type=float, default=WATRAD,
            help='Minimum distance between solute and solvent atoms')
    parser.add_argument('--solvent-selection', default=WATSEL,
            help='Selection of solvent atoms used in solvent-solute distance check')
    parser.add_argument('--solvent-dist', type=float, default=WATCON,
            help='Minimum distance between periodic contacts of solvent atoms')
    parser.add_argument('-v', '--verbose', default=False, action='store_true',
            help='Add chatter to stdout')
    return parser

