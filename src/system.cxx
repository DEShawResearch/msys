#include "system.hxx"
#include "term_table.hxx"
#include <sstream>
#include <stack>
#include <stdexcept>
#include <boost/foreach.hpp>

using namespace desres::msys;

IdList System::_empty;

System::System() 
: _atomprops(ParamTable::create()), _bondprops(ParamTable::create()) {
    //printf("hello structure %p\n", (void *)this);
}

System::~System() {
    //printf("goodbye structure %p\n", (void *)this);
}

Id System::addAtom(Id residue) { 
    Id id = _atoms.size();
    atom_t atm;
    atm.residue = residue;
    atm.gid = id;
    _residueatoms.at(residue).push_back(id);
    _atoms.push_back(atm);
    _atomprops->addParam();
    while (_bondindex.size() < _atoms.size()) {
        _bondindex.push_back(IdList());
    }
    return id;
}

template <typename T>
static IdList get_ids( const T& list, const IdSet& dead ) {
    IdList ids(list.size()-dead.size());
    if (!dead.size()) for (Id i=0; i<ids.size(); i++) {
        ids[i] = i;
    } else {
        Id j=0;
        for (Id i=0; i<list.size(); i++) {
            if (!dead.count(i)) {
                ids[j++]=i;
            }
        }
    }
    return ids;
}

/* remove id from ids.  Assumes id appears at most once in ids */
static void find_and_remove(IdList& ids, Id id) {
    IdList::iterator r = std::remove(ids.begin(), ids.end(), id);
    if (r!=ids.end()) ids.pop_back();
}

IdList System::atoms() const {
    return get_ids(_atoms, _deadatoms);
}
IdList System::bonds() const {
    return get_ids(_bonds, _deadbonds);
}
IdList System::residues() const {
    return get_ids(_residues, _deadresidues);
}
IdList System::chains() const {
    return get_ids(_chains, _deadchains);
}

Id System::findBond( Id i, Id j) const {
    Id n = _bondindex.size();
    if (i>=n || j>=n) return BadId;
    if (i>j) std::swap(i,j);
    const IdList& s = _bondindex[i].size() < _bondindex[j].size() ?
                      _bondindex[i] : _bondindex[j];
    for (IdList::const_iterator b=s.begin(), e=s.end(); b!=e; ++b) {
        const bond_t& bond = _bonds[*b];
        if (bond.i==i && bond.j==j) return *b;
    }
    return BadId;
}

Id System::addBond(Id i, Id j) { 
    if (i>j) std::swap(i,j);
    Id id = findBond(i,j);
    if (!bad(id)) return id;

    if (i==j || j>=_bondindex.size()) {
        std::stringstream ss;
        ss << "addBond: invalid atom ids " << i << ", " << j;
        throw std::runtime_error(ss.str());
    }

    id = _bonds.size();
    _bonds.push_back(bond_t(i,j));
    _bondindex[i].push_back(id);
    _bondindex[j].push_back(id);
    _bondprops->addParam();
    return id;
}

Id System::addResidue(Id chain) {
    Id id = _residues.size();
    residue_t v;
    v.chain = chain;
    _chainresidues.at(chain).push_back(id);
    _residues.push_back(v);
    _residueatoms.push_back(IdList());
    return id;
}

Id System::addChain() {
    Id id = _chains.size();
    chain_t v;
    _chains.push_back(v);
    _chainresidues.push_back(IdList());
    return id;
}

void System::delBond(Id id) {
    const bond_t& b = _bonds.at(id);
    _deadbonds.insert(id);
    find_and_remove(_bondindex[b.i], id);
    find_and_remove(_bondindex[b.j], id);
}

void System::delAtom(Id id) {
    if (id>=_atoms.size()) return;
    IdList del = bondsForAtom(id);
    for (IdList::const_iterator i=del.begin(); i!=del.end(); ++i) {
        delBond(*i);
    }
    _deadatoms.insert(id);
    find_and_remove(_residueatoms.at(_atoms[id].residue), id);
    for (TableMap::iterator t=_tables.begin(); t!=_tables.end(); ++t) {
        t->second->delTermsWithAtom(id);
    }
}

void System::setResidue(Id atm, Id res) {
    Id oldres = _atoms.at(atm).residue;
    if (oldres == res) return;
    /* remove from previous residue */
    find_and_remove(_residueatoms.at(oldres), atm);
    _residueatoms.at(res).push_back(atm);
    _atoms.at(atm).residue = res;
}

void System::delResidue(Id id) {
    /* nothing to do if invalid residue */
    if (id>=_residues.size()) return;

    /* nothing to do if already deleted */
    if (!_deadresidues.insert(id).second) return;

    /* remove from parent chain */
    find_and_remove(_chainresidues.at(_residues[id].chain), id);

    /* remove child atoms.  Clear the index first to avoid O(N) lookups */
    IdList ids;
    ids.swap(_residueatoms.at(id));
    for (IdList::const_iterator iter=ids.begin(); iter!=ids.end(); ++iter) {
        delAtom(*iter);
    }
}

void System::delChain(Id id) {
    /* nothing to do if invalid chain */
    if (id>=_chains.size()) return;

    /* nothing to do if already deleted */
    if (!_deadchains.insert(id).second) return;

    /* remove child residues.  Clear the index first to avoid O(N) lookups */
    IdList ids;
    ids.swap(_chainresidues.at(id));
    for (IdList::const_iterator iter=ids.begin(); iter!=ids.end(); ++iter) {
        delResidue(*iter);
    }
}

Id System::updateFragids(MultiIdList* fragments) {

    /* Create local storage for all atoms (even deleted)
     * this simplifies and speeds up the code below. */
    IdList assignments(_atoms.size(),BadId);

    std::stack<Id> S;
    Id fragid=0;
    for (Id idx=0, n=_atoms.size(); idx<n; idx++) {
        if (!bad(assignments[idx])) continue;   /* already assigned */
        if (_deadatoms.count(idx)) continue;  /* ignore removed atoms */
        S.push(idx);
        assignments[idx] = fragid;
        do {
            Id aid=S.top();
            S.pop();
            IdList& bonds = _bondindex[aid];
            for (IdList::iterator j=bonds.begin(), e=bonds.end(); j!=e; ++j) {
                Id other = _bonds[*j].other(aid);
                if(bad(assignments[other])){
                   assignments[other]=fragid;
                   /* Only add this atom if its non-terminal */
                   if(_bondindex[other].size() >1) S.push(other);
                }
            }
        } while (S.size());
        ++fragid;
    }
    
    /* assign computed fragids to atoms and fragment list */
    if(fragments){
        fragments->clear();
        fragments->resize(fragid);
        for (Id aid=0, n=_atoms.size(); aid<n; aid++) {
            Id fid=assignments[aid];
            _atoms[aid].fragid = fid;
            if(bad(fid)) continue;   /* deleted atom */
            (*fragments)[fid].push_back(aid);
        }
    }else{
        for (Id aid=0, n=_atoms.size(); aid<n; aid++) {
            _atoms[aid].fragid = assignments[aid];
        }
    }
    return fragid;
}

void System::reassignGids() {
    Id gid=0;
    for (Id c=0; c<_chains.size(); c++) {
        if (_deadchains.count(c)) continue;
        BOOST_FOREACH(Id r, _chainresidues[c]) {
            if (_deadresidues.count(r)) continue;
            /* Give pseudos a gid adjacent to their parents.  On the first 
             * pass through the atom list, consider only pseudos.  Make
             * a map from parent atom to pseudo.  On the second pass,
             * when a parent atom is encountered, number its pseudos
             * before proceeding to the next atom.  Pseudos in the residue
             * that aren't bonded to any real atom within the residue
             * will be given gids at the end of the residue's range. */
            std::map<Id,IdList> pseudos;
            std::vector<Id> lone_pseudos;
            BOOST_FOREACH(Id a, _residueatoms[r]) {
                if (_deadatoms.count(a)) continue;
                if (_atoms[a].atomic_number==0) {
                    Id parent = BadId;
                    BOOST_FOREACH(Id b, _bondindex[a]) {
                        Id other = _bonds[b].other(a);
                        if (_atoms[other].residue==r &&
                            _atoms[other].atomic_number>0) {
                            parent = other;
                            break;
                        }
                    }
                    if (bad(parent)) {
                        lone_pseudos.push_back(a);
                    } else {
                        pseudos[parent].push_back(a);
                    }
                }
            }
            BOOST_FOREACH(Id a, _residueatoms[r]) {
                if (_deadatoms.count(a)) continue;
                if (_atoms[a].atomic_number==0) continue;
                _atoms[a].gid=gid++;
                std::map<Id,IdList>::const_iterator plist=pseudos.find(a);
                if (plist!=pseudos.end()) {
                    BOOST_FOREACH(Id p, plist->second) {
                        _atoms[p].gid=gid++;
                    }
                }
            }
            BOOST_FOREACH(Id p, lone_pseudos) {
                _atoms[p].gid=gid++;
            }
        }
    }
}

TermTablePtr System::addTable(const String& name, Id natoms,
                                 ParamTablePtr params) {
    TermTablePtr terms = table(name);
    std::stringstream ss;
    if (!terms) {
        terms.reset(new TermTable(this->shared_from_this(), natoms, params));
        _tables[name]=terms;
    } else if (terms->atomCount()!=natoms) {
        ss << "Could not add table '" << name << "' with " << natoms 
           << " atoms because a table with the same name and " 
           << terms->atomCount() << " already exists.";
        throw std::runtime_error(ss.str());
    } else if ((params && (terms->params()!= params))) {
        ss << "Could not add table '" << name 
           << "' with explicit param table because a table with the same name"
           << " already exists.";
        throw std::runtime_error(ss.str());
    }
    return terms;
}


std::vector<String> System::tableNames() const {
    std::vector<String> s;
    for (TableMap::const_iterator i=_tables.begin(), e=_tables.end(); i!=e; ++i) {
        s.push_back(i->first);
    }
    return s;
}

String System::tableName(boost::shared_ptr<TermTable const> table) const {
    if (table->system() != shared_from_this()) {
        throw std::runtime_error(
                "tableName: table does not belong to this System");
    }
    for (TableMap::const_iterator i=_tables.begin(), e=_tables.end(); i!=e; ++i) {
        if (i->second==table) return i->first;
    }
    /* um, that's bad. */
    assert(false);
    return "";
}

TermTablePtr System::table(const String& name) const {
    TableMap::const_iterator i=_tables.find(name);
    if (i==_tables.end()) return TermTablePtr();
    return i->second;
}

void System::delTable(const String& name) {
    _tables.erase(name);
}

void System::removeTable(TermTablePtr terms) {
    for (TableMap::iterator i=_tables.begin(), e=_tables.end(); i!=e; ++i) {
        if (i->second==terms) {
            _tables.erase(i);
            return;
        }
    }
}

void System::addExtra(const String& name, ParamTablePtr ptr) {
    _extras[name] = ptr;
}

std::vector<String> System::extraNames() const {
    std::vector<String> s;
    for (ExtraMap::const_iterator i=_extras.begin(); i!=_extras.end(); ++i) {
        s.push_back(i->first);
    }
    return s;
}

ParamTablePtr System::extra(const String& name) const {
    ExtraMap::const_iterator i=_extras.find(name);
    if (i==_extras.end()) return ParamTablePtr();
    return i->second;
}

void System::delExtra(const String& name) {
    _extras.erase(name);
}

void System::removeExtra(ParamTablePtr extra) {
    for (ExtraMap::iterator i=_extras.begin(); i!=_extras.end(); ++i) {
        if (i->second==extra) {
            _extras.erase(i);
            return;
        }
    }
}

boost::shared_ptr<System> System::create() {
    return boost::shared_ptr<System>(new System);
}

IdList System::bondedAtoms(Id id) const {
    const IdList& bonds = _bondindex.at(id);
    unsigned i,n = bonds.size();
    IdList ids(n, BadId);
    for (i=0; i<n; i++) ids[i] = bond(bonds[i]).other(id);
    return ids;
}

Id System::maxAtomId() const {
    return _atoms.size();
}
Id System::maxBondId() const {
    return _bonds.size();
}
Id System::maxResidueId() const {
    return _residues.size();
}
Id System::maxChainId() const {
    return _chains.size();
}

Id System::atomPropCount() const {
    return _atomprops->propCount();
}

String System::atomPropName(Id i) const {
    return _atomprops->propName(i);
}

ValueType System::atomPropType(Id i) const {
    return _atomprops->propType(i);
}

Id System::atomPropIndex(String const& name) const {
    return _atomprops->propIndex(name);
}

Id System::addAtomProp(String const& name, ValueType type) {
    return _atomprops->addProp(name,type);
}
void System::delAtomProp(Id index) {
    _atomprops->delProp(index);
}

ValueRef System::atomPropValue(Id term, Id index) {
    return _atomprops->value(term, index);
}

ValueRef System::atomPropValue(Id term, String const& name) {
    return _atomprops->value(term, atomPropIndex(name));
}

Id System::bondPropCount() const {
    return _bondprops->propCount();
}

String System::bondPropName(Id i) const {
    return _bondprops->propName(i);
}

ValueType System::bondPropType(Id i) const {
    return _bondprops->propType(i);
}

Id System::bondPropIndex(String const& name) const {
    return _bondprops->propIndex(name);
}

Id System::addBondProp(String const& name, ValueType type) {
    return _bondprops->addProp(name,type);
}
void System::delBondProp(Id index) {
    _bondprops->delProp(index);
}

ValueRef System::bondPropValue(Id term, Id index) {
    return _bondprops->value(term, index);
}

ValueRef System::bondPropValue(Id term, String const& name) {
    return _bondprops->value(term, bondPropIndex(name));
}

