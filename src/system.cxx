#include "system.hxx"
#include "append.hxx"
#include <msys/version.hxx>
#include <sstream>
#include <stack>
#include <stdexcept>
#include <stdio.h>
#include <ctype.h>
#ifdef DESRES_OS_Darwin
#ifndef _DARWIN_C_SOURCE
#define _DARWIN_C_SOURCE
#endif
#endif
#ifndef _MSC_VER
#include <sys/time.h>
#endif

using namespace desres::msys;

IdList System::_empty;

System::System() 
: _atomprops(ParamTable::create()), _bondprops(ParamTable::create()) {
}

System::~System() {
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
IdList System::cts() const {
    return get_ids(_cts, _deadcts);
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

Id System::addChain(Id ct) {
    if (bad(ct)) {
        if (ctCount()) {
            ct = cts().at(0);
        } else {
            ct = addCt();
        }
    }
    Id id = _chains.size();
    chain_t v;
    v.ct = ct;
    _ctchains.at(ct).push_back(id);
    _chains.push_back(v);
    _chainresidues.push_back(IdList());
    return id;
}

Id System::addCt() {
    Id id = _cts.size();
    _cts.push_back(component_t());
    _ctchains.push_back(IdList());
    return id;
}

component_t::component_t() {
    _kv = ParamTable::create();
    _kv->addParam();
    _kv->addProp("msys_name", StringType);
}

String component_t::name() const { return _kv->value(0,0); }
void component_t::setName(String const& s) { _kv->value(0,0) = s; }

std::vector<String> component_t::keys() const {
    std::vector<String> k;
    for (Id i=1; i<_kv->propCount(); i++) {
        k.push_back(_kv->propName(i));
    }
    return k;
}

void component_t::del(String const& key) { _kv->delProp(_kv->propIndex(key)); }

Id component_t::add(String const& key, ValueType type) {
    return _kv->addProp(key,type);
}

ValueType component_t::type(String const& key) const { 
    return _kv->propType(_kv->propIndex(key));
}

bool component_t::has(String const& key) const {
    return !bad(_kv->propIndex(key));
}

ValueRef component_t::value(String const& key) {
    return _kv->value(0,key);
}
ValueRef component_t::value(Id key) {
    return _kv->value(0,key);
}



component_t& component_t::operator=(component_t const& c) {
    _kv = ParamTable::create();
    AppendParams(_kv, c._kv, c._kv->params());
    return *this;
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

void System::setChain(Id res, Id chn) {
    Id oldchn = _residues.at(res).chain;
    if (oldchn == chn) return;
    /* remove from previous chain */
    find_and_remove(_chainresidues.at(oldchn), res);
    _chainresidues.at(chn).push_back(res);
    _residues.at(res).chain = chn;
}

void System::setCt(Id chn, Id ct) {
    Id oldct = _chains.at(chn).ct;
    if (oldct == ct) return;
    /* remove from previous ct */
    find_and_remove(_ctchains.at(oldct), chn);
    _ctchains.at(ct).push_back(chn);
    _chains.at(chn).ct = ct;
}

IdList System::atomsForCt(Id ct) const {
    IdList ids;
    for (Id const& chn : chainsForCt(ct)) {
        for (Id const& res : residuesForChain(chn)) {
            IdList const& atms = atomsForResidue(res);
            ids.insert(ids.end(), atms.begin(), atms.end());
        }
    }
    // ensure we don't get atoms in a different order just because we
    // accessed them by ct.
    std::sort(ids.begin(), ids.end());
    return ids;
}

Id System::atomCountForCt(Id ct) const {
    Id n=0;
    for (Id const& chn : chainsForCt(ct)) {
        for (Id const& res : residuesForChain(chn)) {
            n += atomCountForResidue(res);
        }
    }
    return n;
}

IdList System::bondsForCt(Id ct) const {
    IdList ids;
    for (Id const& chn : chainsForCt(ct)) {
        for (Id res : residuesForChain(chn)) {
            IdList const& atms = atomsForResidue(res);
            for (Id atm : atms) {
                IdList const& bonds = bondsForAtom(atm);
                for (Id bnd : bonds) {
                    if (bond(bnd).i == atm) ids.push_back(bnd);
                }
            }
        }
    }
    // ensure we don't get bonds in a different order just because we
    // accessed them by ct.
    std::sort(ids.begin(), ids.end());
    return ids;
}

Id System::bondCountForCt(Id ct) const {
    Id n=0;
    for (Id const& chn : chainsForCt(ct)) {
        for (Id res : residuesForChain(chn)) {
            for (Id atm : atomsForResidue(res)) {
                n += bondCountForAtom(atm);
            }
        }
    }
    return n/2;
}

IdList System::residuesForCt(Id ct) const {
    IdList ids;
    for (Id const& chn : chainsForCt(ct)) {
        IdList const& reslist = residuesForChain(chn);
        ids.insert(ids.end(), reslist.begin(), reslist.end());
    }
    return ids;
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

    /* remove from parent ct */
    find_and_remove(_ctchains.at(_chains[id].ct), id);

    /* remove child residues.  Clear the index first to avoid O(N) lookups */
    IdList ids;
    ids.swap(_chainresidues.at(id));
    for (IdList::const_iterator iter=ids.begin(); iter!=ids.end(); ++iter) {
        delResidue(*iter);
    }
}

void System::delCt(Id id) {
    /* nothing to do if invalid ct */
    if (id>=_cts.size()) return;

    /* nothing to do if already deleted */
    if (!_deadcts.insert(id).second) return;

    /* remove child chains .  Clear the index first to avoid O(N) lookups */
    IdList ids;
    ids.swap(_ctchains.at(id));
    for (IdList::const_iterator iter=ids.begin(); iter!=ids.end(); ++iter) {
        delChain(*iter);
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

IdList System::orderedIds() const {
    IdList ids;
    for (Id c=0; c<_chains.size(); c++) {
        if (_deadchains.count(c)) continue;
        for (Id r : _chainresidues[c]) {
            if (_deadresidues.count(r)) continue;
            /* Give pseudos an id adjacent to their parents.  On the first 
             * pass through the atom list, consider only pseudos.  Make
             * a map from parent atom to pseudo.  On the second pass,
             * when a parent atom is encountered, number its pseudos
             * before proceeding to the next atom.  Pseudos in the residue
             * that aren't bonded to any real atom within the residue
             * will be given gids at the end of the residue's range. */
            std::map<Id,IdList> pseudos;
            std::vector<Id> lone_pseudos;
            for (Id a : _residueatoms[r]) {
                if (_deadatoms.count(a)) continue;
                if (_atoms[a].atomic_number==0) {
                    Id parent = BadId;
                    for (Id b : _bondindex[a]) {
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
            for (Id a : _residueatoms[r]) {
                if (_deadatoms.count(a)) continue;
                if (_atoms[a].atomic_number==0) continue;
                ids.push_back(a);
                //_atoms[a].gid=gid++;
                std::map<Id,IdList>::const_iterator plist=pseudos.find(a);
                if (plist!=pseudos.end()) {
                    for (Id p : plist->second) {
                        //_atoms[p].gid=gid++;
                        ids.push_back(p);
                    }
                }
            }
            for (Id p : lone_pseudos) {
                //_atoms[p].gid=gid++;
                ids.push_back(p);
            }
        }
    }
    return ids;
}

TermTablePtr System::addTable(const String& name, Id natoms,
                                 ParamTablePtr params) {
    if (name.empty()) {
        MSYS_FAIL("Table names must have at least one character");
    }
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

void System::renameTable(String const& oldname, String const& newname) {
    TableMap::iterator i=_tables.find(oldname);
    if (i==_tables.end()) {
        std::stringstream ss;
        ss << "renameTable: No such table '" << oldname << "'";
        throw std::runtime_error(ss.str());
    }
    if (_tables.find(newname)!=_tables.end()) {
        std::stringstream ss;
        ss << "renameTable: table named '" << oldname << "' already exists";
        throw std::runtime_error(ss.str());
    }
    TermTablePtr t = i->second;
    _tables.erase(i);
    _tables[newname] = t;
}


String System::tableName(std::shared_ptr<TermTable const> table) const {
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
    TableMap::iterator it = _tables.find(name);
    if (it==_tables.end()) return;
    TermTablePtr t = it->second;
    _tables.erase(it);
    t->destroy();
}

void System::removeTable(TermTablePtr terms) {
    for (TableMap::iterator i=_tables.begin(), e=_tables.end(); i!=e; ++i) {
        if (i->second==terms) {
            delTable(i->first);
            return;
        }
    }
}

void System::coalesceTables() {
    for (TableMap::iterator i=_tables.begin(), e=_tables.end(); i!=e; ++i) {
        i->second->coalesce();
    }
}

void System::addAuxTable(String const& name, ParamTablePtr aux) {
    _auxtables[name] = aux;
}

std::vector<String> System::auxTableNames() const {
    std::vector<String> s;
    AuxTableMap::const_iterator i;
    for (i=_auxtables.begin(); i!=_auxtables.end(); ++i) {
        s.push_back(i->first);
    }
    return s;
}

ParamTablePtr System::auxTable(String const& name) const {
    AuxTableMap::const_iterator i=_auxtables.find(name);
    if (i==_auxtables.end()) return ParamTablePtr();
    return i->second;
}

void System::delAuxTable(String const& name) {
    _auxtables.erase(name);
}

void System::removeAuxTable(ParamTablePtr aux) {
    for (AuxTableMap::iterator i=_auxtables.begin(); i!=_auxtables.end(); ++i) {
        if (i->second==aux) {
            _auxtables.erase(i);
            return;
        }
    }
}

SystemPtr System::create() {
    SystemPtr sys(new System);
    return sys;
}

IdList System::bondedAtoms(Id id) const {
    const IdList& bonds = _bondindex.at(id);
    unsigned i,n = bonds.size();
    IdList ids(n, BadId);
    for (i=0; i<n; i++) ids[i] = bond(bonds[i]).other(id);
    return ids;
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
    /* disallow properties that would conflict with existing columns */
    static const char* badprops[] = {
        "id", "anum", "name", "x", "y", "z", "vx", "vy", "vz",
        "resname", "resid", "chain", "segid", "mass", "charge"
    };
    for (unsigned i=0; i<sizeof(badprops)/sizeof(badprops[0]); i++) {
#ifdef _MSC_VER
        if (!stricmp(name.c_str(), badprops[i])) {
#else
        if (!strcasecmp(name.c_str(), badprops[i])) {
#endif
            MSYS_FAIL("Could not add atom property '" << name << "'\n"
            << "because it would conflict with an existing atom, residue, or chain property");
        }
    }
    return _atomprops->addProp(name,type);
}
void System::delAtomProp(Id index) {
    _atomprops->delProp(index);
}

ValueRef System::atomPropValue(Id term, Id index) {
    return _atomprops->value(term, index);
}

ValueRef System::atomPropValue(Id term, String const& name) {
    Id col = atomPropIndex(name);
    if (bad(col)) MSYS_FAIL("No such atom property '" << name << "'");
    return atomPropValue(term, col);
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
    Id col = bondPropIndex(name);
    if (bad(col)) MSYS_FAIL("No such bond property '" << name << "'");
    return bondPropValue(term, col);
}

double desres::msys::now() {
//FIXME: windows msys, use boost::chrono
#ifndef _MSC_VER
  struct timeval tm;
  struct timezone tz;

  gettimeofday(&tm, &tz);
  return((double)(tm.tv_sec) + (double)(tm.tv_usec)/1000000.0);
#else
  return 0.0;
#endif
}

#define MERGE_NONBONDED_INFO(field) do { \
    if (field.size() && nb.field.size() && field != nb.field) { \
        MSYS_FAIL("incompatible " << field << ": '" << field \
        << "' != '" << nb.field << "'"); \
    } \
    if (field.empty()) field = nb.field; \
} while (0)
void NonbondedInfo::merge(NonbondedInfo const& nb) {
    MERGE_NONBONDED_INFO(vdw_funct);
    MERGE_NONBONDED_INFO(vdw_rule);
    MERGE_NONBONDED_INFO(es_funct);
}
#undef MERGE_NONBONDED_INFO

void GlobalCell::merge(GlobalCell const& gc) {
    double zero[9] = {0,0,0, 0,0,0, 0,0,0};
    const size_t sz = sizeof(zero);
    double* oldvec = (*this)[0];
    const double* newvec = gc[0];
    if (!memcmp(oldvec,zero,sz)) memcpy(oldvec, newvec,sz);
}

int desres::msys::abi_version() {
    return MSYS_ABI_VERSION;
}

