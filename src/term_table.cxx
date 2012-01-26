#include "term_table.hxx"
#include "system.hxx"
#include <stdexcept>
#include <sstream>
#include <iostream>

using namespace desres::msys;

void TermTable::incref(Id p) {
    if (bad(p)) return;
    Id sz = _paramrefs.size();
    _paramrefs.resize(std::max(sz, p+1), 0);
    ++_paramrefs[p];
}

void TermTable::decref(Id p) {
    if (bad(p)) return;
    Id sz = _paramrefs.size();
    _paramrefs.resize(std::max(sz, p+1), BadId);
    --_paramrefs[p];
}

Id TermTable::paramRefs(Id param) const {
    if (param>=_paramrefs.size()) return 0;
    return _paramrefs[param];
}

TermTable::TermTable( SystemPtr system, Id natoms, ParamTablePtr ptr ) 
: _system(system), _natoms(natoms), _props(ParamTable::create()) {
    _params = ptr ? ptr : ParamTable::create();
    _params->incref();
}

TermTable::~TermTable() {
    _params->decref();
}

String TermTable::name() const {
    return system()->tableName(shared_from_this());
}

IdList TermTable::terms() const {
    IdList ids(termCount());
    unsigned i,j=0,n=termCount() + _deadterms.size();
    for (i=0; i<n; i++) {
        if (!_deadterms.count(i)) ids[j++] = i;
    }
    return ids;
}

bool TermTable::alchemical() const {
    return _paramB.size()>0;
}

Id TermTable::termCount() const { 
    return _terms.size()/(1+_natoms) - _deadterms.size(); 
}

Id TermTable::addTerm(const IdList& atoms, Id param) {
    if (atoms.size() != _natoms) {
        throw std::runtime_error("addTerm: incorrect atom count");
    }
    System const& sys = *system();
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
    incref(param);
    _props->addParam();
    return id;
}

void TermTable::delTerm(Id id) {
    if (!hasTerm(id)) return;
    decref(param(id));
    decref(paramB(id));
    _deadterms.insert(id);
    _paramB.erase(id);
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
    decref(this->param(term));
            _terms.at((1+term)*(1+_natoms)-1) = param;
    incref(param);
}

Id TermTable::paramB(Id term) const {
    if (!hasTerm(term)) throw std::runtime_error("Invalid term");
    ParamMap::const_iterator i=_paramB.find(term);
    if (i==_paramB.end()) return BadId;
    return i->second;
}
void TermTable::setParamB(Id term, Id param) {
    if (!hasTerm(term)) throw std::runtime_error("Invalid term");
    ParamMap::const_iterator i=_paramB.find(term);

    Id oldparam = i==_paramB.end() ? BadId : i->second;

    if (bad(param)) {
        decref(oldparam);
        _paramB.erase(term);
    } else if (!_params->hasParam(param)) {
        throw std::runtime_error("Invalid param");
    } else {
        decref(oldparam);
        _paramB[term]=param;
        incref(param);
    }
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

ValueRef TermTable::propValueB(Id term, Id index) {
    Id row = this->paramB(term);
    return _params->value(row, index);
}

ValueRef TermTable::propValueB(Id term, String const& name) {
    Id row = this->paramB(term);
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
    Id i,n = maxTermId();
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
        /* we clear the list, rather than erase the entire map, so that
         * we can reuse the set of unique params for paramB */
        ids.clear();
    }
    /* map paramB for each term to a distinct set, using the params from
     * the first pass when possible */
    for (i=0; i<n; i++) {
        if (!hasTerm(i)) continue;
        Id p = paramB(i);
        if (bad(p)) continue;
        map[p].push_back(i);
    }
    for (ParamMap::iterator iter=map.begin(); iter!=map.end(); ++iter) {
        Id p = iter->first;
        IdList& ids = iter->second;
        for (Id i=0; i<ids.size(); i++) setParamB(ids[i], p);
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

