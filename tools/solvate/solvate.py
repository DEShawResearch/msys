
import os, msys

WAT=os.path.join(os.getenv('MSYS_PREFIX'), 'share', 'solvate', 'h2o.dms')
WATRAD = 2.4

def Solvate(mol, thickness=0, minmax=None,
            xpad=None, ypad=None, zpad=None, 
            chain='WT', verbose=False,
            cubic=False):
    '''
    Build a solvent shell around the given solute.  If no solute is provided,
    the solute is treated as a point at the origin.  thickness specifies the
    amount of solute added along each coordinate axis; this can be overridden
    by xpad, ypad, zpad.  The extent of the solvent can also be given 
    explicitly with minmax.  If cubic is true, the box size will be expanded
    to the size of the longest dimension.
    '''

    wat=msys.LoadDMS(WAT)
    watsize=[wat.cell[i][i] for i in range(3)]

    print "finding bounding box..."
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

    if verbose: print "replicating %d x %d x %d" % (nx,ny,nz)
    return
    wat=cPickle.load(file(WAT))
    watchain=mol.addChain(chain)
    # construct the template water
    res=watchain.appendResidue()
    res.name='H2O'
    o, h1, h2 = res.appendAtoms(3)
    o.name="O"
    h1.name="H1"
    h2.name="H2"
    o.num=8
    h1.num=1
    h2.num=1
    mol.addBond(o,h1)
    mol.addBond(o,h2)
    # replicate the template water
    nwat = len(wat)/3 * nx * ny * nz - 1
    if verbose: print "appending %d residues..." % nwat
    watchain.appendResidues( nwat, res )
    # map the coordinates to the waters
    iter=watchain.atoms.__iter__()
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                shift=[xmin+i*watsize[0],ymin+j*watsize[1],zmin+k*watsize[2]]
                coords = wat + shift
                for pos in coords:
                  a = iter.next()
                  a.pos = pos


    view = mol.select('chain %s and same residue as (x<%f or y<%f or z<%f or x>%f or y>%f or z>%f or within %f of (not chain %s))' % (
        chain,xmin,ymin,zmin,xmax,ymax,zmax,WATRAD, chain))
    if verbose: 
      print "removing %d atoms, %d residues" % (
            len(view.atoms), len(view.residues) )
    mol.removeAtoms(view.atomList())
    mol.prune()
    if verbose: print "updating global cell to (%g %g %g)" % (
        xmax-xmin, ymax-ymin, zmax-zmin)

    # update the cell
    mol.cell[0]=geom.Vec3(xmax-xmin,0,0)
    mol.cell[1]=geom.Vec3(0,ymax-ymin,0)
    mol.cell[2]=geom.Vec3(0,0,zmax-zmin)

