#include "system.hxx"
#include "sssr.hxx"
#include "term_table.hxx"
#include <sstream>
#include <stack>
#include <queue>
#include <stdexcept>
#include <boost/foreach.hpp>
#include <boost/algorithm/string.hpp> /* for boost::trim */
#include <stdio.h>
#ifdef DESRES_OS_Darwin
#ifndef _DARWIN_C_SOURCE
#define _DARWIN_C_SOURCE
#endif
#endif
#include <sys/time.h>

using namespace desres::msys;

IdList System::_empty;

namespace {
  struct atomsel_macro {
      const char * name;
      const char * text;
  };

  atomsel_macro builtin_macros[] = {
      {"at","resname ADE A THY T"},
      {"acidic","resname ASP GLU"},
      {"cyclic","resname HIS PHE PRO TRP TYR"},
      {"acyclic","protein and not cyclic"},
      {"aliphatic","resname ALA GLY ILE LEU VAL"},
      {"alpha","protein and name CA"},
      {"amino","protein"},
      {"aromatic","resname HIS PHE TRP TYR"},
      {"basic","resname ARG HIS LYS HSP"},
      {"bonded","numbonds > 0"},
      {"buried"," resname ALA LEU VAL ILE PHE CYS MET TRP"},
      {"cg","resname CYT C GUA G"},
      {"charged","basic or acidic"},
      {"hetero","not (protein or nucleic)"},

    // APD's hydrophobic residue list, from Branden and Tooze (pp6-7).
      {"hydrophobic","resname ALA LEU VAL ILE PRO PHE MET TRP"},
    
      {"small","resname ALA GLY SER"},
      {"medium","resname VAL THR ASP ASN PRO CYS ASX PCA HYP"},
      {"large","protein and not (small or medium)"},
      {"neutral","resname VAL PHE GLN TYR HIS CYS MET TRP ASX GLX PCA HYP"},
      {"polar","protein and not hydrophobic"},
      {"purine","resname ADE A GUA G"},
      {"pyrimidine","resname CYT C THY T URA U"},
      {"surface","protein and not buried"},
      {"lipid","resname DLPE DMPC DPPC GPC LPPC PALM PC PGCL POPC POPE"},
      {"lipids","lipid"},
      {"ion","resname AL BA CA Ca CAL CD CES CLA CL Cl CO CS CU Cu CU1 CUA HG IN IOD K MG MN3 MO3 MO4 MO5 MO6 NA Na NAW OC7 PB POT PT RB SOD TB TL WO4 YB ZN ZN1 ZN2"},
      {"ions","ion"},
      {"sugar","resname AGLC"},
      {"solvent","not (protein or sugar or nucleic or lipid)"},
      /* for carbon, nitrogen, oxygen, and sulfur (and hydrogen), VMD
       * uses a name-based regex; e.g. 'name "N.*"' for nitrogen.  This
       * is just silly, and gets things like Na and Cl wrong.  We refuse
       * to reproduce this buggy behavior and intead look to the atomic
       * number. */
      {"carbon","atomicnumber 6"},
      {"nitrogen","atomicnumber 7"},
      {"oxygen","atomicnumber 8"},
      {"sulfur","atomicnumber 16"},
      {"noh","not hydrogen"},
      {"heme","resname HEM HEME"}
  };
}

void System::initSelectionMacros() {
    const unsigned n=sizeof(builtin_macros)/sizeof(builtin_macros[0]);
    for (unsigned i=0; i<n; i++) {
        addSelectionMacro(builtin_macros[i].name, builtin_macros[i].text);
    }
}

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
    BOOST_FOREACH(Id const& chn, chainsForCt(ct)) {
        BOOST_FOREACH(Id const& res, residuesForChain(chn)) {
            IdList const& atms = atomsForResidue(res);
            ids.insert(ids.end(), atms.begin(), atms.end());
        }
    }
    return ids;
}

Id System::atomCountForCt(Id ct) const {
    Id n=0;
    BOOST_FOREACH(Id const& chn, chainsForCt(ct)) {
        BOOST_FOREACH(Id const& res, residuesForChain(chn)) {
            n += atomCountForResidue(res);
        }
    }
    return n;
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
        BOOST_FOREACH(Id r, _chainresidues[c]) {
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
                ids.push_back(a);
                //_atoms[a].gid=gid++;
                std::map<Id,IdList>::const_iterator plist=pseudos.find(a);
                if (plist!=pseudos.end()) {
                    BOOST_FOREACH(Id p, plist->second) {
                        //_atoms[p].gid=gid++;
                        ids.push_back(p);
                    }
                }
            }
            BOOST_FOREACH(Id p, lone_pseudos) {
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
    sys->initSelectionMacros();
    return sys;
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
Id System::maxCtId() const {
    return _cts.size();
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
        if (!strcasecmp(name.c_str(), badprops[i])) {
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

void System::addSelectionMacro(std::string const& macro,
                               std::string const& definition) {
    _macros[macro]=definition;
}

std::string const& System::selectionMacroDefinition(std::string const& m) const {
    MacroMap::const_iterator it=_macros.find(m);
    if (it==_macros.end()) {
        static const std::string _empty;
        return _empty;
    }
    return it->second;
}

Id System::selectionMacroCount() const {
    return _macros.size();
}

std::vector<std::string> System::selectionMacros() const {
    std::vector<std::string> v;
    for (MacroMap::const_iterator it=_macros.begin(); it!=_macros.end(); ++it) {
        v.push_back(it->first);
    }
    return v;
}

void System::delSelectionMacro(std::string const& macro) {
    _macros.erase(macro);
}

void System::clearSelectionMacros() {
    _macros.clear();
}

void System::copySelectionMacros(System const& m) {
    _macros = m._macros;
}

void SystemImporter::initialize(IdList const& atoms) {
    resmap.clear();
    chnmap.clear();
    IdList reslist, chnlist;
    BOOST_FOREACH(Id atm, atoms) {
        Id res = sys->atom(atm).residue;
        Id chn = sys->residue(res).chain;
        reslist.push_back(res);
        chnlist.push_back(chn);
    }
    sort_unique(reslist);
    sort_unique(chnlist);

    BOOST_FOREACH(Id chn, chnlist) {
        chain_t const& chain = sys->chain(chn);
        chnmap[ChnKey(chain.ct, chain.name, chain.segid)] = chn;
        BOOST_FOREACH(Id res, sys->residuesForChain(chn)) {
            if (std::binary_search(reslist.begin(), reslist.end(), res)) {
                residue_t const& residue = sys->residue(res);
                resmap[ResKey(chn, residue.resid, residue.name, residue.insertion)] = res;
            }
        }
    }
}

Id SystemImporter::addAtom(std::string chain, std::string segid,
                           int resnum, std::string resname,
                           std::string aname,
                           std::string insertion,
                           Id ct) {

    boost::trim(chain);
    boost::trim(segid);
    boost::trim(resname);
    boost::trim(insertion);
    boost::trim(aname);

    if (bad(ct)) MSYS_FAIL("Got ct=BadId");

    /* start a new ct if necessary */
    while (!sys->hasCt(ct)) {
        sys->addCt();
    }

    /* start a new chain if necessary */
    std::pair<ChnMap::iterator,bool> cp;
    cp = chnmap.insert(std::make_pair(ChnKey(ct,chain,segid),chnid));
    if (cp.second) {
        /* new chain/segid pair found, so start a new chain */
        chnid = sys->addChain(ct);
        sys->chain(chnid).name = chain;
        sys->chain(chnid).segid = segid;
        resid = BadId;
        cp.first->second = chnid;
    } else {
        /* use existing chain */
        chnid = cp.first->second;
    }

    /* start a new residue if necessary */
    std::pair<ResMap::iterator,bool> p;
    p = resmap.insert(std::make_pair(ResKey(chnid,resnum,resname,insertion), 
                resid));
    if (p.second) {
        /* new resname/resnum in this chain, so start a new residue. */
        resid = sys->addResidue(chnid);
        sys->residue(resid).resid = resnum;
        sys->residue(resid).name = resname;
        sys->residue(resid).insertion = insertion;
        p.first->second = resid;
    } else {
        /* use existing residue */
        resid = p.first->second;
    }

    /* add the atom */
    Id atm = sys->addAtom(resid);
    sys->atom(atm).name = aname;
    return atm;
}

Id System::glueCount() const {
    return _glue.size();
}

std::vector<glue_t> System::gluePairs() const {
    std::vector<glue_t> v(_glue.size());
    std::copy(_glue.begin(), _glue.end(), v.begin());
    return v;
}

bool System::hasGluePair(Id p0, Id p1) const {
    if (p0>p1) std::swap(p0,p1);
    return _glue.count(glue_t(p0,p1));
}

bool System::addGluePair(Id p0, Id p1) {
    if (p0>p1) std::swap(p0,p1);
    if (!hasAtom(p0)) MSYS_FAIL("Bad atom id: " << p0);
    if (!hasAtom(p1)) MSYS_FAIL("Bad atom id: " << p1);
    if (p0==p1)       MSYS_FAIL("Identical glue atoms with id " << p0);
    return _glue.insert(glue_t(p0,p1)).second;
}

bool System::delGluePair(Id p0, Id p1) {
    if (p0>p1) std::swap(p0,p1);
    return _glue.erase(glue_t(p0,p1));
}


namespace {

    bool is_water( System* sys, Id res ) {
        Id O(BadId), H1(BadId), H2(BadId);
        IdList const& atoms = sys->atomsForResidue(res);
        for (Id i=0; i<atoms.size(); i++) {
            Id id=atoms[i];
            const atom_t& b = sys->atom(id);
            if (b.atomic_number==8) {
                if (bad(O)) O=id;
                else {
                    O=BadId;
                    break;
                }
            } else if (b.atomic_number==1) {
                if      (bad(H1)) H1=id;
                else if (bad(H2)) H2=id;
                else {
                    O=BadId;
                    break;
                }
            }
        }
        if (bad(O) || bad(H1) || bad(H2)) return false;
        return
            sys->bondCountForAtom(O)==2 &&
            sys->bondCountForAtom(H1)==1 &&
            sys->bondCountForAtom(H2)==1 &&
            sys->bond(sys->bondsForAtom(H1)[0]).other(H1)==O &&
            sys->bond(sys->bondsForAtom(H2)[0]).other(H2)==O;
    }

    bool has_water_residue_name( const std::string& resname ) {
        static const char * names[] = {
            "H2O", "HH0", "OHH", "HOH", "OH2", "SOL", "WAT", "TIP",
            "TIP2", "TIP3", "TIP4", "SPC"
        };
        unsigned i,n = sizeof(names)/sizeof(names[0]);
        for (i=0; i<n; i++) {
            if (resname==names[i]) return true;
        }
        return false;
    }

    void find_sidechain(System* mol, Id res, Id ca) {
        /* pick a c-beta atom, or settle for a hydrogen */
        Id cb=BadId;
        BOOST_FOREACH(Id nbr, mol->bondedAtoms(ca)) {
            if (mol->atom(nbr).type==AtomProBack) continue;
            if (bad(cb)) {
                cb=nbr;
            } else {
                /* already have a cb candidate.  Pick the better one. */
                if (mol->atom(nbr).atomic_number==1 || 
                    toupper(mol->atom(nbr).name[0]=='H')) continue;
                cb=nbr;
            }
        }
        if (bad(cb)) return;
        /* starting from cb, recursively add all atoms in the residue
         * which are bonded to the cb but are not backbone. */
        std::queue<Id> q;
        q.push(cb);
        while (!q.empty()) {
            Id atm = q.front();
            mol->atom(atm).type = AtomProSide;
            BOOST_FOREACH(Id nbr, mol->bondedAtoms(atm)) {
                atom_t const& nbratm = mol->atom(nbr);
                if (nbratm.type==AtomOther && nbratm.residue==res) {
                    q.push(nbr);
                }
            }
            q.pop();
        }
    }

    void analyze_residue(System* sys, Id res) {

        static const char * protypes[] = { "CA", "C", "O", "N" };
        /* NOTE: VMD 1.7 _tries_ to include O2 in its backbone selection, but
         * because of a bug in the implementation, it doesn't include it.  
         * It's fixed in 1.9.1. and in vmd/1.9.0-17.  */
        static const char * proterms[] = {
            "OT1", "OT2", "OXT", "O1", "O2"
        };
        static const char * nuctypes[] = {
            "P", "O1P", "O2P", "OP1", "OP2", "C3*", "C3'", "O3*", "O3'",
            "C4*", "C4'", "C5*", "C5'", "O5*", "O5'"
        };
        static const char * nucterms[] = {
            "H5T", "H3T"
        };
        typedef std::map<std::string,AtomType> NameMap;
        NameMap types, terms;
        if (types.empty()) {
            for (unsigned i=0; i<sizeof(protypes)/sizeof(protypes[0]); i++) {
                types[protypes[i]]=AtomProBack;
            }
            for (unsigned i=0; i<sizeof(nuctypes)/sizeof(nuctypes[0]); i++) {
                types[nuctypes[i]]=AtomNucBack;
            }
            for (unsigned i=0; i<sizeof(proterms)/sizeof(proterms[0]); i++) {
                terms[proterms[i]]=AtomProBack;
            }
            for (unsigned i=0; i<sizeof(nucterms)/sizeof(nucterms[0]); i++) {
                terms[nucterms[i]]=AtomNucBack;
            }
        }

        /* clear structure */
        IdList const& atoms = sys->atomsForResidue(res);
        for (Id i=0; i<atoms.size(); i++) sys->atom(atoms[i]).type=AtomOther;
        sys->residue(res).type = ResidueOther;

        /* check for water */
        if (has_water_residue_name(sys->residue(res).name) || 
                is_water(sys,res)) {
            sys->residue(res).type = ResidueWater;
            return;
        }

        /* need at least four atoms to determine protein or nucleic */
        if (atoms.size()<4) return;

        int npro=0, nnuc=0;
        std::set<std::string> names;

        Id ca_atm = BadId;
        Id c_atm = BadId;
        Id n_atm = BadId;
        for (Id i=0; i<atoms.size(); i++) {
            Id id = atoms[i];
            const atom_t& atm = sys->atom(id);
            const std::string& aname = atm.name;
            if (!names.insert(aname).second) continue;
            /* check for nucleic or protein backbone */
            NameMap::const_iterator iter=types.find(aname);
            AtomType atype=AtomOther;
            if (iter!=types.end()) {
                atype=iter->second;
                if (atype==AtomProBack) {
                    if (aname=="CA") ca_atm = id;
                    else if (aname=="C") c_atm = id;
                    else if (aname=="N") n_atm = id;
                }
            } else {
                /* try terminal names */
                iter=terms.find(aname);
                if (iter!=terms.end()) {
                    /* must be bonded to atom of the same type */
                    IdList const& bonds = sys->bondsForAtom(id);
                    for (Id j=0; j<bonds.size(); j++) {
                        Id o = sys->bond(bonds[j]).other(id);
                        AtomType otype = sys->atom(o).type;
                        if (otype==iter->second) {
                            atype=otype;
                            break;
                        }
                    }
                }
            }
            sys->atom(id).type = atype;
            if (atype==AtomProBack) ++npro;
            if (atype==AtomNucBack) ++nnuc;
        }
        ResidueType rtype=ResidueOther;
        if (npro>=4 && 
            ca_atm!=BadId && c_atm!=BadId && n_atm != BadId &&
            !bad(sys->findBond(ca_atm,n_atm)) &&
            !bad(sys->findBond(ca_atm,c_atm)) &&
             bad(sys->findBond(c_atm,n_atm))) {
            rtype=ResidueProtein;
        } else if (nnuc>=4) {
            rtype=ResidueNucleic;
        } else for (Id i=0; i<atoms.size(); i++) {
            sys->atom(atoms[i]).type = AtomOther;
        }
        if (rtype==ResidueProtein) {
            find_sidechain(sys, res, ca_atm);
        }
        sys->residue(res).type=rtype;
    }
}

void System::analyze() {
    updateFragids();
    for (Id res=0; res<maxResidueId(); res++) {
        if (hasResidue(res)) analyze_residue(this,res);
    }
}

double desres::msys::now() {
  struct timeval tm;
  struct timezone tz;

  gettimeofday(&tm, &tz);
  return((double)(tm.tv_sec) + (double)(tm.tv_usec)/1000000.0);
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
    const Vec3 zero;
    if (A==zero && B==zero && C==zero) {
        *this = gc;
    } else if (A!=gc.A || B!=gc.B || C != gc.C) {
        MSYS_FAIL("Incompatible global cells");
    }
}
