"""
dms-thermalize input.dms output.dms [ options ]

Assign Boltzmann-sampled velocities to the atoms.  Atoms with zero mass
will get zero velocity.

"""
import msys
import random, math
import numpy
import sys, os

# the Boltzmann constant in J/K
BOLTZMANN = 1.3806503e-23
AVOGADRO = 6.0221415e23


def apply(mol, T, seed=None):
    """assign random velocities sampled from a Boltzmann distribution
    of temperature T.
    """
    g = random
    if seed is not None:
        g = g.Random(seed)

    # convert to amu * (A/ps)^2
    kT = BOLTZMANN * T * AVOGADRO / 10.0
    for a in mol.atoms:
        m = a.mass
        if m == 0:
            t = 0.0
        else:
            t = math.sqrt(kT / m)
        a.vx = t * g.gauss(0, 1)
        a.vy = t * g.gauss(0, 1)
        a.vz = t * g.gauss(0, 1)


def remove_drift(mol):
    """Remove center of mass motion.  Returns the original center of mass
    velocity.   Zero out the velocity of pseudoparticles.
    """
    vel = mol.getVelocities()
    mass = [a.mass for a in mol.atoms]
    avg = numpy.average(vel, axis=0, weights=mass)
    vel -= avg
    vel[mol.selectIds("atomicnumber 0")] = 0
    mol.setVelocities(vel)
    return avg


def main():
    import optparse

    parser = optparse.OptionParser(__doc__)
    parser.add_option(
        "-t",
        "--temperature",
        type="float",
        default=300,
        help="thermalization temperature in Kelvin",
    )
    parser.add_option(
        "-s",
        "--seed",
        default="1",
        help="Random seed, or 'random' to get a random random seed",
    )
    parser.add_option(
        "--remove-drift",
        action="store_true",
        default=True,
        help="Remove center of mass motion from assigned velocities",
    )
    parser.add_option(
        "--no-remove-drift",
        action="store_false",
        dest="remove_drift",
        help="Do not remove center of mass motion from assigned velocities",
    )

    opts, args = parser.parse_args()
    if len(args) != 2:
        parser.error("incorrect number of arguments")

    print("Loading input file <%s>" % args[0])
    mol = msys.LoadDMS(args[0])

    print("Thermalizing to %sK using seed %s" % (opts.temperature, opts.seed))

    seed = opts.seed
    if seed == "random":
        seed = None
    else:
        seed = int(seed)

    apply(mol, opts.temperature, seed)

    if opts.remove_drift:
        print("Removing center of mass motion")
        remove_drift(mol)

    print("Writing DMS file <%s>" % args[1])
    msys.SaveDMS(mol, args[1])
