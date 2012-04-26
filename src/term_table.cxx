#include "term_table.hxx"
#include "system.hxx"
#include <stdexcept>
#include <sstream>
#include <iostream>

using namespace desres::msys;

TermTable::TermTable( SystemPtr system, Id natoms, ParamTablePtr ptr ) 
: _system(system), _natoms(natoms), _props(ParamTable::create()),
  category(NO_CATEGORY) {
    _params = ptr ? ptr : ParamTable::create();
}

void TermTable::destroy() {
    SystemPtr sys = system();
    if (!sys) return;
    _system.reset();
    sys->removeTable(shared_from_this());
    for (unsigned i=0; i<maxTermId(); i++) {
        if (_deadterms.count(i)) continue;
        _params->decref(param(i));
    }
}

String TermTable::name() const {
    SystemPtr sys = system();
    if (!sys) MSYS_FAIL("Table has been destroyed");
    return sys->tableName(shared_from_this());
}

void TermTable::rename(String const& newname) {
    SystemPtr sys = system();
    if (!sys) MSYS_FAIL("Table has been destroyed");
    return sys->renameTable(name(), newname);
}

IdList TermTable::terms() const {
    IdList ids(termCount());
    unsigned i,j=0,n=termCount() + _deadterms.size();
    for (i=0; i<n; i++) {
        if (!_deadterms.count(i)) ids[j++] = i;
    }
    return ids;
}

Id TermTable::termCount() const { 
    return _terms.size()/(1+_natoms) - _deadterms.size(); 
}

Id TermTable::addTerm(const IdList& atoms, Id param) {
    if (atoms.size() != _natoms) {
        throw std::runtime_error("addTerm: incorrect atom count");
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
    Id id=_terms.size()/(1+_natoms);
    _terms.insert(_terms.end(), atoms.begin(), atoms.end());
    _terms.push_back(param);
    _params->incref(param);
    _props->addParam();
    return id;
}

void TermTable::delTerm(Id id) {
    if (!hasTerm(id)) return;
    _params->decref(param(id));
    _deadterms.insert(id);
}

void TermTable::delTermsWithAtom(Id atm) {
    Id stride = 1+_natoms;
    Id i,n = _terms.size()/stride;
    for (i=0; i<n; i++) {
        if (_deadterms.count(i)) continue;
        TermList::const_iterator t=_terms.begin() + i*stride;
        TermList::const_iterator e=t+_natoms;
        if (std::find(t, e, atm)!=e) delTerm(i);
    }
}

Id TermTable::param(Id term) const { 
    if (!hasTerm(term)) throw std::runtime_error("Invalid term");
    return  _terms.at((1+term)*(1+_natoms)-1);
}

void TermTable::setParam(Id term, Id param) {
    if (!hasTerm(term)) throw std::runtime_error("Invalid term");
    if (!(bad(param) || _params->hasParam(param))) {
        MSYS_FAIL("Invalid param " << param << " for table " << name());
    }
    _params->decref(this->param(term));
    _terms.at((1+term)*(1+_natoms)-1) = param;
    _params->incref(param);
}

bool TermTable::hasTerm(Id term) const {
    return term<_terms.size() && !_deadterms.count(term);
}

Id TermTable::maxTermId() const {
    return _terms.size()/(1+_natoms);
}

IdList TermTable::atoms(Id i) const { 
    TermList::const_iterator b=_terms.begin()+i*(1+_natoms);
    return IdList(b, b+_natoms);
}

Id TermTable::atom(Id term, Id index) const {
    TermList::const_iterator b=_terms.begin()+term*(1+_natoms);
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
    return _props->value(term, index);
}

ValueRef TermTable::termPropValue(Id term, String const& name) {
    return _props->value(term, termPropIndex(name));
}

ValueRef TermTable::propValue(Id term, Id index) {
    Id row = this->param(term);
    return _params->value(row, index);
}

ValueRef TermTable::propValue(Id term, String const& name) {
    Id row = this->param(term);
    return _params->value(row, name);
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
        for (Id i=0; i<ids.size(); i++) setParam(ids[i], p);
    }
}

static const char* category_names[] = {
    "none",
    "bond",
    "constraint",
    "virtual",
    "polar",
    "nonbonded",
    "exclusion"
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

