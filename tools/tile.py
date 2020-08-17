"""
Tile a system system by nx, ny, nz
"""
import msys
import numpy
import itertools


def parse_args():
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("inputs", nargs="+", type=msys.Load, help="input file(s)")
    parser.add_argument("-o", "--output", required=True, help="output file")
    parser.add_argument("--nx", default=1, type=int, help="tiles in the x direction")
    parser.add_argument("--ny", default=1, type=int, help="tiles in the y direction")
    parser.add_argument("--nz", default=1, type=int, help="tiles in the z direction")
    parser.add_argument(
        "--allow-unequal-cells",
        action="store_true",
        default=False,
        help="allow input files to have different cell shapes",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    nx, ny, nz = args.nx, args.ny, args.nz
    if nx < 1 or ny < 1 or nz < 1:
        print("Tile dimensions must be positive")
        return 1
    n = nx * ny * nz

    mols = args.inputs
    if n % len(mols) != 0:
        print(
            "Number of input files (%d) is incommensurate with number of replicas (%d)"
            % (len(mols), n)
        )
        return 1

    cell = mols[0].getCell()
    for i, mol in enumerate(mols):
        pos = mol.positions
        posrange = pos.max(0) - pos.min(0)
        for celldim, extent in zip(numpy.diag(cell), posrange):
            if celldim * 1.5 < extent:
                print(
                    "Cannot tile system %d: box size appears to be much smaller than"
                    " spatial extent of particles"
                    % (i + 1)
                )
                return 1
        icell = mols[i].getCell().tolist()
        if cell.tolist() != icell:
            print(
                "The global cell of the %d'th molecule differs from the first" % (i + 1)
            )
            print(cell.tolist(), "!=", icell)
            if not args.allow_unequal_cells:
                return 1

    cycle = itertools.cycle(mols)

    print("Tiling %d x %d x %d: %d total copies" % (nx, ny, nz, n))
    xshift = -0.5 * (nx - 1) * cell[0]
    yshift = -0.5 * (ny - 1) * cell[1]
    zshift = -0.5 * (nz - 1) * cell[2]

    out = msys.CreateSystem()

    for i in range(nx):
        xdelta = xshift + i * cell[0]
        for j in range(ny):
            ydelta = yshift + j * cell[1]
            for k in range(nz):
                zdelta = zshift + k * cell[2]
                delta = xdelta + ydelta + zdelta
                mol = next(cycle)
                mol.translate(delta)
                out._ptr.append(mol._ptr, msys._msys.BadId)
                mol.translate(-delta)

    # set up the unit cell
    out.setCell(numpy.dot(numpy.diag((nx, ny, nz)), cell))

    # copy nonbonded info
    out.nonbonded_info = mol.nonbonded_info

    out.coalesceTables()
    out = out.clone()

    print("Saving '%s'" % args.output)
    msys.Save(out, args.output)

    print("Done")
