from . import _msys


class experimental:
    @staticmethod
    def make_rules(rules, flavor):
        if flavor != "amber":
            raise ValueError("unsupported flavor '%s'; must be 'amber'" % flavor)
        rules.vdw_func = "vdw_12_6"
        rules.vdw_rule = "arithmetic/geometric"
        rules.exclusions = 4
        rules.es_scale = [0.0, 0.0, 0.8333]
        rules.lj_scale = [0.0, 0.0, 0.5]
        return rules

    @staticmethod
    def build_pairs_and_exclusions(mol, flavor="amber"):
        """Build all pair_12_6_es and exclusion terms for the given system

        Arguments:
            mol (msys.System): chemical system
            flavor (str): 'amber' for Amber-style exclusion and 1-4 scaling

        Notes:
            Currently does not support systems with pseudoatoms (virtuals or drudes).
            This is checked.

        """
        if mol.selectIds("mass <= 0"):
            raise ValueError("Systems with pseudos not yet supported")

        tuples = _msys.Tuples()
        for frag in mol.updateFragids():
            tuples.build(mol._ptr, [a.id for a in frag])

        ff = _msys.Forcefield()
        experimental.make_rules(ff.rules, flavor)
        ff.build_component(_msys.Component.exclusions, mol._ptr, tuples)
