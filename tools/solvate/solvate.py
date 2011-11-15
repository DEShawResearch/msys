
import os, msys

WAT=os.path.join(os.getenv('MSYS_PREFIX'), 'share', 'solvate', 'h2o.dms')
WATRAD = 2.4

def Solvate(mol, watbox=None, dims=None, center=None,
        chain='WT', verbose=False):
    '''
    Build and return a new system consisting of mol plus solvent.

    If no box dimensions are provided, the box size along each axis will
    be equal to the greater of that of the input structure and the input
    water box.

    Return the solvated system; no modifications are made to the input system.
    '''

    mol=mol.clone()
    if watbox is None: watbox=WAT

    if verbose:
        print "Loading water box from '%s'" % watbox
    wat=msys.LoadDMS(watbox)
    watsize=[wat.cell[i][i] for i in range(3)]
    molsize=[mol.cell[i][i] for i in range(3)]

    # find a chain name for the waters that doesn't overlap with the 
    # input structure.
    chains = set(c.name for c in mol.chains)
    watchain = 'X'
    while watchain in chains:
        watchain += 'X'
    for c in wat.chains:
        c.name = watchain

    if verbose:
        print "Dims specified as '%s'" % dims
    if dims is None:
        dims = [max(x) for x in zip(watsize, molsize)]
    elif len(dims)==1:
        dims=[float(dims[0])]*3
    elif len(dims)!=3:
        raise ValueError, "Dims must be given as list of one or three values"
    else:
        dims=[float(x) for x in dims]
    
    if center is None:
        center=[0,0,0]
    elif len(center)!=3:
        raise ValueError, "Center must be given as list of three values"
    else:
        center=[float(x) for x in center]

    if center==[0,0,0]:
        WITHIN='within'
    else:
        WITHIN='pbwithin'


    # update the cell
    mol.cell[0][:]=[dims[0],0,0]
    mol.cell[1][:]=[0,dims[1],0]
    mol.cell[2][:]=[0,0,dims[2]]

    # extract where to put the water
    xmin = center[0]-0.5*dims[0]
    ymin = center[1]-0.5*dims[1]
    zmin = center[2]-0.5*dims[2]
    xmax = center[0]+0.5*dims[0]
    ymax = center[1]+0.5*dims[1]
    zmax = center[2]+0.5*dims[2]
    nx = int(dims[0]/watsize[0]) + 1
    ny = int(dims[1]/watsize[1]) + 1
    nz = int(dims[2]/watsize[2]) + 1

    xshift = -0.5 * (nx-1)*watsize[0]
    yshift = -0.5 * (ny-1)*watsize[1]
    zshift = -0.5 * (nz-1)*watsize[2]

    # replicate the template water box
    if verbose: print "replicating %d x %d x %d" % (nx,ny,nz)
    for i in range(nx):
        xdelta = xshift + i*watsize[0]
        for j in range(ny):
            ydelta = yshift + j*watsize[1]
            for k in range(nz):
                zdelta = zshift + k*watsize[2]
                newatoms = mol.append(wat)
                for a in newatoms:
                    a.x += xdelta
                    a.y += ydelta
                    a.z += zdelta

    if verbose: print "removing overlaps"

    toofar='x<%s or y<%s or z<%s or x>%s or y>%s or z>%s' % (
            xmin,ymin,zmin,xmax,ymax,zmax)

    toonear='%s %s of not chain %s' % (WITHIN, WATRAD, watchain)
    mol=mol.clone(
            'not same residue as (chain %s and oxygen and (%s or %s))' % (
                watchain, toofar, toonear))

    # assign the water chain name and water resids
    watres = 1
    for c in mol.chains:
        if c.name == watchain:
            c.name = chain
            for r in c.residues:
                r.resid = watres
                watres += 1

    mol.reassignGids()

    if verbose: print "updating global cell to (%g %g %g)" % tuple(dims)

    return mol

