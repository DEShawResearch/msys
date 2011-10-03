#include "system.hxx"
#include "term_table.hxx"
#include <sstream>
#include <stack>
#include <stdexcept>

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

bool System::findBond( Id i, Id j, Id * id ) const {
    Id n = _bondindex.size();
    if (i>=n || j>=n) return false;
    if (i>j) std::swap(i,j);
    const IdList& s = _bondindex[i].size() < _bondindex[j].size() ?
                      _bondindex[i] : _bondindex[j];
    for (IdList::const_iterator b=s.begin(), e=s.end(); b!=e; ++b) {
        const bond_t& bond = _bonds[*b];
        if (bond.i==i && bond.j==j) {
            if (id) *id = *b;
            return true;
        }
    }
    return false;
}

Id System::addBond(Id i, Id j) { 
    Id id = _bonds.size();
    if (i>j) std::swap(i,j);
    if (findBond(i,j,&id)) return id;
    if (i==j || j>=_bondindex.size()) return BadId;
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

Id System::update_fragids(MultiIdList* fragments) {

    /* Create local storage for all atoms (even deleted)
     * this simplifies and speeds up the code below. */
    IdList assignments(_atoms.size(),BadId);

    std::stack<Id> S;
    Id fragid=0;
    for (Id aid=0, n=_atoms.size(); aid<n; aid++) {
        if (!bad(assignments[aid])) continue;   /* already assigned */
        if (_deadatoms.count(aid)) continue;  /* ignore removed atoms */
        S.push(aid);
        assignments[aid] = fragid;
        do {
            aid=S.top();
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
    } else if ((params && (terms->paramTable()!= params))) {
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

Id System::renumberGids(Id starting_gid) {
    IdList ids=atoms();
    for (IdList::const_iterator i=ids.begin(), e=ids.end(); i!=e; ++i) {
        atom(*i).gid = starting_gid++;
    }
    return starting_gid;
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

