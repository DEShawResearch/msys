import msys


def guess_type(val):
    try:
        int(val)
    except ValueError:
        pass
    else:
        return int

    try:
        float(val)
    except ValueError:
        pass
    else:
        return float

    return str


def dms_set_atoms(atoms, attr, val):
    if not len(atoms):
        return
    if hasattr(msys.Atom, attr):
        val = type(getattr(atoms[0], attr))(val)
        for atom in atoms:
            setattr(atom, attr, val)
    else:
        t = guess_type(val)
        val = t(val)
        atoms[0].system.addAtomProp(attr, t)
        for atom in atoms:
            atom[attr] = val


def dms_set_residues(atoms, attr, val):
    if not len(atoms):
        return
    residues = set()
    val = type(getattr(atoms[0].residue, attr))(val)
    for atom in atoms:
        res = atom.residue
        if res not in residues:
            setattr(res, attr, val)
            residues.add(res)


def dms_set_chains(atoms, attr, val):
    if not len(atoms):
        return
    chains = set()
    val = type(getattr(atoms[0].residue.chain, attr))(val)
    for atom in atoms:
        chn = atom.residue.chain
        if chn not in chains:
            setattr(chn, attr, val)
            chains.add(chn)


def dms_set_box(mol, attr, val):
    val = float(val)
    box = mol.cell
    if attr == "d":
        for i in range(3):
            box[i][i] = val
    else:
        try:
            i = {"x": 0, "y": 1, "z": 2}[attr]
        except KeyError:
            raise ValueError("Box attribute must be x, y, z or d, got '%s'" % attr)
        box[i][i] = val


def dms_set_term(atoms, table, attr, val):
    if attr in table.term_props:
        val = table.termPropType(attr)(val)
    elif attr in table.params.props:
        val = table.params.propType(attr)(val)
    else:
        raise ValueError("No property '%s' in table '%s'" % (attr, table.name))

    atoms = set(atoms)
    for t in table.terms:
        if atoms.issuperset(t.atoms):
            t[attr] = val
    table.coalesce()


def Update(mol, atoms, key, val):
    """update the system by setting properties corresponding to key to the
    value val.  key can take the following forms:
    1) foo          -- same as atom.foo
    2) atom.foo     -- atom property foo
    3) residue.foo  -- residue property foo
    4) chain.foo    -- chain property foo
    5) table.foo    -- property foo of all terms
    6) box.foo      -- where foo is x, y, z, or d (sets all three)
    """
    ndots = key.count(".")
    if ndots == 0:
        key = ".".join(("atom", key))
    elif ndots > 1:
        raise ValueError("key '%s' can have no more than one '.'" % key)
    target, prop = key.split(".")
    if target == "atom":
        dms_set_atoms(atoms, prop, val)
    elif target == "residue":
        dms_set_residues(atoms, prop, val)
    elif target == "chain":
        dms_set_chains(atoms, prop, val)
    elif target == "box":
        dms_set_box(mol, prop, val)
    elif target in mol.table_names:
        dms_set_term(atoms, mol.table(target), prop, val)
    else:
        raise ValueError("no table named '%s'" % target)
