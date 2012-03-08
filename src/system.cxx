#include "system.hxx"
#include "term_table.hxx"
#include "atomsel/vmd.hxx"
#include <sstream>
#include <stack>
#include <stdexcept>
#include <boost/foreach.hpp>

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
        //addSelectionMacro(builtin_macros[i].name, builtin_macros[i].text);
        /* The addSelectionMacro method checks whether the macro is valid
         * before adding, which takes extra time.  We don't need to check
         * the built-in macros. */
        _macros[builtin_macros[i].name] = builtin_macros[i].text;
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

void System::setChain(Id res, Id chn) {
    Id oldchn = _residues.at(res).chain;
    if (oldchn == chn) return;
    /* remove from previous chain */
    find_and_remove(_chainresidues.at(oldchn), res);
    _chainresidues.at(chn).push_back(res);
    _residues.at(res).chain = chn;
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

void System::addSelectionMacro(std::string const& macro,
                               std::string const& definition) {

    /* prevent recursive macro definitions */
    std::string olddef = selectionMacroDefinition(macro);
    delSelectionMacro(macro);
    try {
        atomsel::vmd::parse(definition, shared_from_this());
    }
    catch (std::runtime_error& e) {
        if (olddef.size()) _macros[macro]=olddef;
        throw;
    }
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
