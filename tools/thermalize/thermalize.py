
import msys
import random, math

# the Boltzmann constant in J/K
BOLTZMANN = 1.3806503e-23
AVOGADRO  = 6.0221415e23

def apply(mol, T):
    ''' assign random velocities sampled from a Boltzmann distribution
    of temperature T.
    '''

    # convert to amu * (A/ps)^2
    kT = BOLTZMANN * T * AVOGADRO / 10.0
    for a in mol.atoms:
        m=a.mass
        if m==0:
            t=0.0
        else:
            t=math.sqrt( kT / m )
        a.vx = t*random.gauss(0,1)
        a.vy = t*random.gauss(0,1)
        a.vz = t*random.gauss(0,1)

