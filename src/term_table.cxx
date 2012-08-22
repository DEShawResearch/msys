#include "term_table.hxx"
#include "system.hxx"
#include "override.hxx"
#include <stdexcept>
#include <sstream>
#include <iostream>
#include <algorithm>

using namespace desres::msys;

TermTable::TermTable( SystemPtr system, Id natoms, ParamTablePtr ptr ) 
: _system(system), _natoms(natoms), _ndead(0), _props(ParamTable::create()),
  _maxIndexId(0), category(NO_CATEGORY) {
    _params = ptr ? ptr : ParamTable::create();
    _overrides = OverrideTable::create(_params);
    if (natoms<1) MSYS_FAIL("TermTable must have at least 1 atom");
}

void TermTable::destroy() {
    SystemPtr sys = system();
    if (!sys) return;
    _system.reset();
    sys->removeTable(shared_from_this());
    for (unsigned i=0; i<maxTermId(); i++) {
        if (!_alive(i)) continue;
        _params->decref(param(i));
    }
    _overrides->clear();
    _terms.clear();
    _index.clear();
}

String TermTable::name() const {
    SystemPtr sys = system();
    if (!sys) return "";
    return sys->tableName(shared_from_this());
}

void TermTable::rename(String const& newname) {
    SystemPtr sys = system();
    if (!sys) MSYS_FAIL("Table has been destroyed");
    return sys->renameTable(name(), newname);
}

IdList TermTable::terms() const {
    IdList ids(termCount());
    Id i,j=0,n=maxTermId();
    for (i=0; i<n; i++) {
        if (_alive(i)) ids[j++]=i;
    }
    return ids;
}

Id TermTable::addTerm(const IdList& atoms, Id param) {
    if (atoms.size() != _natoms) {
        MSYS_FAIL("incorrect atom count for TermTable " << name());
    }
    SystemPtr s = system();
    if (!s) MSYS_FAIL("Table has been destroyed");
    System const& sys = *s;
    for (IdList::const_iterator atm=atoms.begin(); atm!=atoms.end(); ++atm) {
        if (!sys.hasAtom(*atm)) {
            std::stringstream ss;
            ss << "addTerm: no such atom " << *atm;
            throw std::runtime_error(ss.str());
        }
    }
    Id id=maxTermId();
    _terms.insert(_terms.end(), atoms.begin(), atoms.end());
    _terms.push_back(param);
    _params->incref(param);
    _props->addParam();
    return id;
}

void TermTable::delTerm(Id id) {
    if (!hasTerm(id)) return;
    _params->decref(param(id));
    if (_index.size()) {
        Id i,n = atomCount();
        for (i=0; i<n; i++) {
            Id atm = atom(id,i);
            IdList& p = _index.at(atm);
            p.resize(std::remove(p.begin(), p.end(), id)-p.begin());
        }
    }
    _terms[id*(1+_natoms)] = BadId; /* mark as dead */
    ++_ndead;
}

void TermTable::delTermsWithAtom(Id atm) {
    IdList ids(1, atm);
    IdList terms = findWithAll(ids);
    for (unsigned i=0; i<terms.size(); i++) delTerm(terms[i]);
}

Id TermTable::param(Id term) const { 
    if (!hasTerm(term)) {
        MSYS_FAIL("Table '" << name() << "' has no term with id " << term);
    }
    return  _terms.at((1+term)*(1+_natoms)-1);
}

void TermTable::setParam(Id term, Id param) {
    if (!hasTerm(term)) {
        MSYS_FAIL("Table '" << name() << "' has no term with id " << term);
    }
    if (!(bad(param) || _params->hasParam(param))) {
        MSYS_FAIL("Invalid param " << param << " for table " << name());
    }
    _params->decref(this->param(term));
    _terms.at((1+term)*(1+_natoms)-1) = param;
    _params->incref(param);
}

IdList TermTable::atoms(Id term) const { 
    if (term>=maxTermId()) {
        MSYS_FAIL("Table '" << name() << "' has no term with id " << term);
    }
    IdList::const_iterator b=_terms.begin()+term*(1+_natoms);
    return IdList(b, b+_natoms);
}

Id TermTable::atom(Id term, Id index) const {
    if (term>=maxTermId()) {
        MSYS_FAIL("Table '" << name() << "' has no term with id " << term);
    }
    if (index>=atomCount()) {
        MSYS_FAIL("Table '" << name() << "' has no atoms for index " << index);
    }
    IdList::const_iterator b=_terms.begin()+term*(1+_natoms);
    return *(b+index);
}

Id TermTable::termPropCount() const {
    return _props->propCount();
}

String TermTable::termPropName(Id i) const {
    return _props->propName(i);
}

ValueType TermTable::termPropType(Id i) const {
    return _props->propType(i);
}

Id TermTable::termPropIndex(String const& name) const {
    return _props->propIndex(name);
}

Id TermTable::addTermProp(String const& name, ValueType type) {
    if (_params->propIndex(name)!=BadId) {
        MSYS_FAIL("TermTable " << this->name() << " already has a param property '" << name << "'");
    }
    return _props->addProp(name,type);
}

void TermTable::delTermProp(Id index) {
    _props->delProp(index);
}

ValueRef TermTable::termPropValue(Id term, Id index) {
    if (!hasTerm(term)) {
        MSYS_FAIL("TermTable " << this->name() << " has no term with id " << term);
    }
    return _props->value(term, index);
}

ValueRef TermTable::termPropValue(Id term, String const& name) {
    Id index = termPropIndex(name);
    if (bad(index)) {
        std::string extra;
        if (!bad(_params->propIndex(name))) {
                extra = ", although a param property with that name does exist.";
        }
        MSYS_FAIL("TermTable " << this->name() << " has no term property '" << name << "'" << extra);
    }
    return termPropValue(term, index);
}

ValueRef TermTable::propValue(Id term, Id index) {
    Id row = this->param(term);
    return _params->value(row, index);
}

ValueRef TermTable::propValue(Id term, String const& name) {
    Id index = _params->propIndex(name);
    if (bad(index)) {
        std::string extra;
        if (!bad(termPropIndex(name))) {
            extra = ", although a term property with that name does exist.";
        }
        MSYS_FAIL("TermTable " << this->name() << " has no param property '" << name << "'" << extra);
    }
    return propValue(term, index);
}

namespace {
    struct ParamComparator {
        ParamTablePtr params;
        
        explicit ParamComparator(ParamTablePtr p) : params(p) {}
        bool operator()(Id const& pi, Id const& pj) {
            return params->compare(pi,pj)<0;
        }
    };
}

void TermTable::coalesce() {
    if (!_params->paramCount()) return;
    typedef std::map<Id,IdList,ParamComparator> ParamMap;
    IdList old2new(_params->paramCount(), BadId);
    ParamComparator comp(_params);
    assert(!comp(0,0));
    ParamMap map(comp);
    
    /* hash all parameters, so that tables that share a param table
     * will wind up using the same params when possible */
    Id i,n = _params->paramCount();
    for (i=0; i<n; i++) map[i];

    n = maxTermId();
    for (i=0; i<n; i++) {
        if (!hasTerm(i)) continue;
        Id p = param(i);
        if (bad(p)) continue;
        map[p].push_back(i);
    }
    for (ParamMap::iterator iter=map.begin(); iter!=map.end(); ++iter) {
        Id p = iter->first;
        IdList& ids = iter->second;
        for (Id i=0; i<ids.size(); i++) {
            Id t = ids[i];
            old2new.at(param(t)) = p;
            setParam(t, p);
        }
    }

    /* coalesce the override table, if it exists.  Reuse the same set of
     * distinct parameters.  If an override points to a parameter that 
     * is unused, it can be safely removed. */
    std::vector<IdPair> L = _overrides->list();
    for (unsigned i=0; i<L.size(); i++) {
        Id p1 = old2new.at(L[i].first);
        Id p2 = old2new.at(L[i].second);
        Id p = _overrides->get(L[i]);
        _overrides->del(L[i]);
        if (bad(p1) || bad(p2)) continue;
        _overrides->set(IdPair(p1,p2), p);
    }
    /* At this point there may be duplicate params in the override
     * params, but I won't bother fixing that now */
}

static const char* category_names[] = {
    "none",
    "bond",
    "constraint",
    "virtual",
    "polar",
    "nonbonded",
    "exclusion",
    "override"
};
static const unsigned ncategories = sizeof(category_names)/sizeof(char *);


Category desres::msys::parse(std::string const& s) {
    for (unsigned i=0; i<ncategories; i++) {
        if (s==category_names[i]) return (Category)i;
    }
    std::stringstream ss;
    ss << "Could not convert '" << s << "' into msys Category";
    throw std::runtime_error(ss.str());
}

std::string desres::msys::print(Category const& c) {
    unsigned n = (unsigned)c;
    if (n>=ncategories) {
        std::stringstream ss;
        ss << "Unrecognized msys category with numerical value " << n;
        throw std::runtime_error(ss.str());
    }
    return category_names[n];
}


void TermTable::update_index() {
    _index.resize(system()->maxAtomId());
    if (_terms.empty()) return;
    const Id natoms=atomCount();
    Id i=_maxIndexId, n=maxTermId();
    for (; i<n; i++) {
        if (!hasTerm(i)) continue;
        for (Id j=0; j<natoms; j++) {
            Id atm = atom(i,j);
            _index.at(atm).push_back(i);
        }
    }
    _maxIndexId = n;
}

IdList TermTable::findWithAll(IdList const& ids) {
    update_index();
    IdList terms;
    for (unsigned i=0; i<ids.size(); i++) {
        IdList const& p = _index.at(ids[i]);
        if (i==0) {
            terms = p;
        } else {
            IdList v(terms.size());
            v.resize(
                    std::set_intersection(
                        p.begin(), p.end(), 
                        terms.begin(), terms.end(),
                        v.begin())-v.begin());
            terms = v;
        }
    }
    return terms;
}

IdList TermTable::findWithAny(IdList const& ids) {
    update_index();
    IdList terms;
    for (unsigned i=0; i<ids.size(); i++) {
        IdList const& p = _index.at(ids[i]);
        if (i==0) {
            terms = p;
        } else {
            IdList v(terms.size()+p.size());
            v.resize(
                    std::set_union(
                        p.begin(), p.end(), 
                        terms.begin(), terms.end(),
                        v.begin())-v.begin());
            terms = v;
        }
    }
    return terms;
}

IdList TermTable::findExact(IdList const& ids) {
    const unsigned natoms = atomCount();
    IdList terms;
    if (ids.size()!=natoms) return terms;

    /* index gets updated by findWithAll() */
    IdList tmp = findWithAll(ids);
    for (unsigned i=0; i<tmp.size(); i++) {
        Id term = tmp[i];
        for (unsigned j=0; j<natoms; j++) {
            if (atom(term,j)!=ids[j]) {
                term = BadId;
                break;
            }
        }
        if (!bad(term)) terms.push_back(term);
    }
    return terms;
}

void TermTable::resetParams(ParamTablePtr params) {
    for (Id i=0; i<maxTermId(); i++) {
        if (_alive(i)) setParam(i,BadId);
    }
    _params = params;
}
