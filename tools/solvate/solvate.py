
import os, msys

WAT=os.path.join(os.getenv('MSYS_PREFIX'), 'share', 'solvate', 'h2o.dms')
WATRAD = 2.4

def Solvate(mol, thickness=0, minmax=None,
            xpad=None, ypad=None, zpad=None, 
            chain='WT', verbose=False,
            cubic=False):
    '''
    Build and return a new system consisting of mol plus solvent.

    If no solute is provided, the solute is treated as a point at
    the origin.  thickness specifies the amount of solute added along
    each coordinate axis; this can be overridden by xpad, ypad, zpad.
    The extent of the solvent can also be given explicitly with minmax.
    If cubic is true, the box size will be expanded to the size of the
    longest dimension.

    Return the solvated system; no modifications are made to the input system.
    '''

    mol=mol.clone()

    wat=msys.LoadDMS(WAT)
    watsize=[wat.cell[i][i] for i in range(3)]

    # find a chain name for the waters that doesn't overlap with the 
    # input structure.
    chains = set(c.name for c in mol.chains)
    watchain = 'X'
    while watchain in chains:
        watchain += 'X'
    for c in wat.chains:
        c.name = watchain

    if verbose: print "finding bounding box..."
    if len(mol.atoms):
        first=mol.atoms[0].pos
        xmin, ymin, zmin = first
        xmax, ymax, zmax = first
        for a in mol.atoms:
            x,y,z = a.pos
            if   x<xmin: xmin=x
            elif x>xmax: xmax=x
            if   y<ymin: ymin=y
            elif y>ymax: ymax=y
            if   z<zmin: zmin=z
            elif z>zmax: zmax=z
    else:
        xmin, ymin, zmin =[wat.cell[i][i] for i in range(3)]
        xmin=ymin=zmin=-t
        xmax=ymax=zmax= t/2-1e15
    smin=[xmin,ymin,zmin]
    smax=[xmax,ymax,zmax]
    
    # The definition of the solvent range is a little tricky.
    # For each direction, if pad is given, then the range is the
    # solute minmax along that axis padded by the specified amount,
    # no matter what thickness and minmax are.  If pad is not given,
    # then the minmax argument is checke first, and if it's also missing,
    # then we use thickness.  
    pad=(xpad, ypad, zpad)
    for i in range(3):
        t=thickness
        if pad[i] is None:
            if minmax:
                smin[i]=minmax[0][i]
                smax[i]=minmax[1][i]
            else:
                smin[i] -= t
                smax[i] += t
        else:
            t=pad[i]
            smin[i] -= t
            smax[i] += t

    if cubic:
        delta=0.5*max(smax[i]-smin[i] for i in range(3))
        center=[0.5*(smax[i]+smin[i]) for i in range(3)]
        smin=[center[i]-delta for i in range(3)]
        smax=[center[i]+delta for i in range(3)]

    # extract where to put the water
    xmin, ymin, zmin = smin
    xmax, ymax, zmax = smax
    dx, dy, dz = [a-b for a,b in zip(smax,smin)]
    nx = int(dx/watsize[0]) + 1
    ny = int(dy/watsize[1]) + 1
    nz = int(dz/watsize[2]) + 1

    # replicate the template water box
    if verbose: print "replicating %d x %d x %d" % (nx,ny,nz)
    for i in range(nx):
        xdelta = xmin + i*watsize[0]
        for j in range(ny):
            ydelta = ymin + j*watsize[1]
            for k in range(nz):
                zdelta = zmin + k*watsize[2]

                newatoms = mol.append(wat)
                for a in newatoms:
                    a.x += xdelta
                    a.y += ydelta
                    a.z += zdelta

    if verbose: print "removing overlaps"
    # FIXME: use the low-level interface for speed
    mol = mol.clone(
        'not (chain %s and same residue as (x<%f or y<%f or z<%f or x>%f or y>%f or z>%f or within %f of (not chain %s)))' % (
            watchain,xmin,ymin,zmin,xmax,ymax,zmax,WATRAD, watchain))

    # assign the water chain name
    for c in mol.chains:
        if c.name == watchain:
            c.name = chain

    mol.reassignGids()

    if verbose: print "updating global cell to (%g %g %g)" % (
        xmax-xmin, ymax-ymin, zmax-zmin)

    # update the cell
    mol.cell[0][:]=[xmax-xmin,0,0]
    mol.cell[1][:]=[0,ymax-ymin,0]
    mol.cell[2][:]=[0,0,zmax-zmin]

    return mol

