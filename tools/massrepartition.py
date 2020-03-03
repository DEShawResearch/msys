'''
   Repartition hydrogen masses
'''
from __future__ import print_function
import msys
import argparse


class HMRActions(argparse.Action):
    def __call__(self, parser, namespace, values, option_string):
        if len(values) != 2: raise RuntimeError(f"I need 2 values, got: {values}")
        sel = values[0]
        mass = None
        try:
            mass = float(values[1])
        except:
            pass
        if mass is None: raise RuntimeError(f"Second argument to {option_string} must be the mass, got: {values[1]}")
        repartition = {"--repartition": True, "-r": True, "--norepartition": False, "-n": False}[option_string]
        namespace.actions.append((repartition, sel, mass))

def parse_args():
    import argparse
    "Parse the arguments"
    examples ='''\n
* To repartition non-water hydrogens:
    dms-hmr input.dms output.dms -r 'not water' 3.0

* To change all hydrogen to deuterium:
    dms-hmr input.dms output.dms -n 'all' 2.014

* To repartition water and just set mass of non-water hydrogens (piana2013)
    dms-hmr input.dms output.dms -r 'water' 4.0 -n 'not water' 4.0

'''
    epilog='''You can use: "{--repartition, -r, --norepartition, -n} SELECTION MASS" any number of times'''
    parser = argparse.ArgumentParser(usage=examples, epilog=epilog)
    parser.add_argument('--norepartition', '-n', action=HMRActions, dest="actions", 
                        nargs=2, metavar=("SELECTION", "MASS"),
                        help='Set mass of hydrogens in SELECTION to MASS but do NOT repartition from the bonded heavy atom')
    parser.add_argument('--repartition', '-r', action=HMRActions, dest="actions",
                        nargs=2, metavar=("SELECTION", "MASS"),
                        help='Set mass of hydrogens in SELECTION to MASS by repartitioning from the bonded heavy atom')
    parser.set_defaults(actions=[])

    parser.add_argument("input_dms", help="Input structure file")
    parser.add_argument("output_dms", help="Output structure file")
    return parser.parse_args()


def adjust_hmasses(mol, sel, hmass, repartition):
    """
    Adjust masses of hydrogens (for example HMR)

    Parameters
    -----------
    mol : msys.System
        Input system
    sel : str
        Selection of hydrogens to adjust masses (e.g. 'all', 'water')
    hmass : float
        New mass for hydrogens
    repartition : bool
        repartition the mass from the heavy atom bonded to the hydrogen

    Returns
    -------
    msys.System
        System with the masses adjusted

    """
    new_mol = mol.clone()
    atomsel = 'hydrogen and (%s)' % sel
    hatoms = new_mol.select(atomsel)
    print(f"  * Setting mass for '{atomsel}' to {hmass} (N={len(hatoms)})")
    warn = False
    for a in hatoms + [b for h in hatoms for b in h.bonded_atoms]:
        dev = abs(a.mass/msys.MassForElement(a.atomic_number) -1)
        if dev > 0.1: 
            warn = True
            break
    if (warn):
       s = f'  * WARNING: It looks like the hydrogen or bonded atom masses in selection "{sel}" were already modified!'
       s += '\n    You may have overlapping selections or are modifing a file that already has hmr applied.'
       s += '\n    Please verify that the final masses are as expected'
       print(s)
    total_mass_before = sum([a.mass for a in new_mol.atoms])
    for a in hatoms:
        old = a.mass
        a.mass = hmass
        if repartition:
            batom = a.bonded_atoms[0]
            batom.mass -= hmass - old
            assert batom.mass > 0, f"Mass for atom {batom.id} ({batom.name}) is < 0 ({batom.mass})"
    total_mass_after = sum([a.mass for a in new_mol.atoms])
    dmass = total_mass_after - total_mass_before
    if repartition:
        assert dmass < 0.001, f"Total mass has changed too much (dm={dmass})"
        print(f"  * Total system mass before and after: {total_mass_after}")
    else:
        print(f"  * Total system mass changed by {dmass}")
    return new_mol


def main():
    args = parse_args()
    input_dms = args.input_dms
    output_dms = args.output_dms

    print(f"Reading system '{input_dms}'")
    mol = msys.Load(input_dms)
    for (repartition, selection, hmass) in args.actions:
        mol = adjust_hmasses(mol, selection, hmass, repartition)
    print(f"Saving system '{output_dms}'")
    msys.Save(mol, output_dms)

