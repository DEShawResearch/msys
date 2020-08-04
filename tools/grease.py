'''
dms-grease input.dms lipid.dms output.dms [ options ]

Tile a lipid bilayer around a solute.

dms-grease builds a new chemical system consisting of the input system
plus a lipid bilayer constructed by tiling *lipid.dms* in the x-y plane.
If the *input.dms* is given as "-", then a pure membrane will be built.

An error will be encountered if only one of *input.dms* and *lipid.dms* 
have forcefield information; this is because Msys refuses to write DMS
files that have only partial information for the nonbonded atom types.
If you don't have forcefield information for one of the input files,
use the *--structure-only* option to ignore the forcefield information
in the one that does.

The global cell of the new system will be orthorhombic and have x and
y dimensions given by the specified size of the membrane, and z dimension
given by the input structure or the lipid membrane template, whichever is
greater.
'''
import sys, os, msys

def Grease(mol, tile, thickness=0.0, xsize=None, ysize=None,
            ctname='grease', verbose=True, 
            square=False):
    '''
    Build and return a new system consisting of mol plus lipid bilayer.
    Tile is the lipid bilayer system to replicate.

    If no solute is provided, the solute is treated as a point at the
    origin.  thickness specifies the amount of solute added along the x
    and y axis.  The dimensions of the bilayer can also be given explicitly
    with dimensions.  If square is true, the box size will be expanded to
    the size of the longest dimension.

    Lipids will be created in a new Ct.  Their chain names will be left the
    same as in the original tile, but the resids will be renumbered to ensure
    uniqueness.

    Return the greased system; no modifications are made to the input system.
    '''

    mol=mol.clone()

    lipsize=[tile.cell[i][i] for i in range(3)]
    lx, ly, lz = lipsize
    if lx<=1 or ly<=1:
        raise ValueError("Lipid tile has missing box dimensions.")

    if xsize is None or ysize is None:
        if verbose: print("finding bounding box...")
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
            if xsize is None: xsize = thickness + xmax - xmin
            if ysize is None: ysize = thickness + ymax - ymin
        else:
            if xsize is None: xsize = thickness
            if xsize is None: ysize = thickness
        
    xsize = float(xsize)
    ysize = float(ysize)

    if square:
        xsize = max(xsize, ysize)
        ysize = xsize

    # extract where to put the lipid
    nx = int(xsize // lipsize[0]) + 1
    ny = int(ysize // lipsize[1]) + 1

    xshift = -0.5 * (nx-1)*lipsize[0]
    yshift = -0.5 * (ny-1)*lipsize[1]
    xmin = -0.5 * xsize
    ymin = -0.5 * ysize
    xmax = -xmin
    ymax = -ymin

    # update the cell
    mol.cell[0][:]=[xsize,0,0]
    mol.cell[1][:]=[0,ysize,0]
    mol.cell[2][:]=[0,0,max(mol.cell[2][2], lipsize[2])]

    # Create a new Ct for the lipids
    ct = mol.addCt()
    ct.name = 'grease'
    ctnum = ct.id + 1
    # replicate the template lipid box
    if verbose: print("replicating %d x %d" % (nx,ny))
    for i in range(nx):
        xdelta = xshift + i*lipsize[0]
        for j in range(ny):
            ydelta = yshift + j*lipsize[1]
            newatoms = ct.append(tile)
            for a in newatoms:
                a.x += xdelta
                a.y += ydelta

    if verbose: print("replicated system contains %d atoms" % mol.natoms)
    if verbose: print("removing overlap with solute")
    headgroup_dist = 2.0
    mol = mol.clone(
        'not (ctnumber %d and same residue as (atomicnumber 8 15 and pbwithin %f of (noh and not ctnumber %d)))' % (
            ctnum, headgroup_dist, ctnum))
    dist = 1.0
    mol = mol.clone(
        'not (ctnumber %d and same residue as (pbwithin %f of (noh and not ctnumber %d)))' % (
            ctnum, dist, ctnum))
    if verbose: print("after removing solute overlap, have %d atoms" % mol.natoms)

    if verbose: print("removing outer lipids and water")
    mol = mol.clone(
        'not (ctnumber %d and same residue as (atomicnumber 8 15 and (abs(x)>%f or abs(y)>%f)))' % (
            ctnum,xmax,ymax))
    if verbose: print("after removing outer lipids, have %d atoms" % mol.natoms)

    # renumber lipid resids
    lipnum = 1
    for c in ct.chains:
            for r in c.residues:
                r.resid = lipnum
                lipnum += 1

    return mol

def main():
    import optparse
    parser = optparse.OptionParser(__doc__)

    parser.add_option('-t', '--thickness', default=0.0, type='float',
        help='Minimum distance from edge of membrane to input structure')
    parser.add_option('-x', '--xsize', default=None,
        help='Size of membrane along x dimension')
    parser.add_option('-y', '--ysize', default=None,
        help='Size of membrane along y dimension')
    parser.add_option(      '--square', default=False, action='store_true',
        help='Ensure the resulting membrane is square')
    parser.add_option('--structure-only', default=False, action='store_true',
        help='Ignore forcefield in input.dms and lipid.dms')
    parser.add_option('-v', '--verbose', default=True, action='store_true',
        help="Be chatty")

    opts, args = parser.parse_args()
    if len(args)!=3:
        parser.error("incorrect number of arguments")

    input, lipid, output = args
    structure_only = opts.structure_only
    del opts.__dict__['structure_only']

    if opts.verbose: print("Loading input file <%s>" % input)
    if input=='-':
        mol=msys.CreateSystem()
    else:
        mol=msys.LoadDMS(input, structure_only=structure_only)

    tile=msys.LoadDMS(lipid, structure_only=structure_only)
    
    mol = Grease(mol, tile, **opts.__dict__)

    if opts.verbose: print("Writing DMS file <%s>" % output)
    msys.SaveDMS(mol,output)

