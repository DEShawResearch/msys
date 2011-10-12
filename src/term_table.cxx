#include "term_table.hxx"
#include <stdexcept>

using namespace desres::msys;

TermTable::TermTable( SystemPtr system, Id natoms, ParamTablePtr ptr ) 
: _system(system), _natoms(natoms), _props(ParamTable::create()) {
    _params = ptr ? ptr : ParamTable::create();
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
    Id id=_terms.size()/(1+_natoms);
    _terms.insert(_terms.end(), atoms.begin(), atoms.end());
    _terms.push_back(param);
    _props->addParam();
    return id;
}

void TermTable::delTerm(Id id) {
    /* don't delete nonexistent term */
    if (id<_terms.size()) _deadterms.insert(id);
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
    return  _terms.at((1+term)*(1+_natoms)-1);
}

void TermTable::setParam(Id term, Id param) {
            _terms.at((1+term)*(1+_natoms)-1) = param;
}

Id TermTable::paramB(Id term) const {
    ParamMap::const_iterator i=_paramB.find(term);
    if (i==_paramB.end()) return BadId;
    return i->second;
}
void TermTable::setParamB(Id term, Id param) {
    if (bad(param)) {
        _paramB.erase(term);
    } else if (!hasTerm(term)) {
        throw std::runtime_error("Invalid term");
    } else if (!_params->hasParam(param)) {
        throw std::runtime_error("Invalid param");
    } else {
        _paramB[term]=param;
    }
}

bool TermTable::hasTerm(Id term) const {
    return term<_terms.size() && !_deadterms.count(term);
}

IdList TermTable::atoms(Id i) const { 
    TermList::const_iterator b=_terms.begin()+i*(1+_natoms);
    return IdList(b, b+_natoms);
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

bool TermTable::findTermProp(String const& name, Id* indexptr) const {
    Id id = _props->propIndex(name);
    if (indexptr) *indexptr = id;
    return !bad(id);
}

Id TermTable::addTermProp(String const& name, ValueType type) {
    return _props->addProp(name,type);
}

ValueRef TermTable::termPropValue(Id term, Id index) {
    return _props->value(term, index);
}

ValueRef TermTable::termPropValue(Id term, String const& name) {
    return _props->value(term, termPropIndex(name));
}

Id TermTable::propCount() const {
    return _params->propCount();
}

Id TermTable::paramCount() const {
    return _params->paramCount();
}

String TermTable::propName(Id i) const {
    return _params->propName(i);
}

ValueType TermTable::propType(Id i) const {
    return _params->propType(i);
}

Id TermTable::propIndex(String const& name) const {
    return _params->propIndex(name);
}

bool TermTable::findProp(String const& name, Id* indexptr) const {
    Id id = _params->propIndex(name);
    if (indexptr) *indexptr = id;
    return !bad(id);
}

Id TermTable::addProp(String const& name, ValueType type) {
    return _params->addProp(name,type);
}

ValueRef TermTable::propValue(Id term, Id index) {
    return _params->value(param(term), index);
}

ValueRef TermTable::propValueB(Id term, Id index) {
    Id paramid = paramB(term);
    if (bad(paramid)) paramid = param(term);
    return _params->value(paramid, index);
}

ValueRef TermTable::propValue(Id term, String const& name) {
    return _params->value(param(term), propIndex(name));
}

ValueRef TermTable::propValueB(Id term, String const& name) {
    Id paramid = paramB(term);
    if (bad(paramid)) paramid = param(term);
    return _params->value(paramid, propIndex(name));
}

