#include "param_table.hxx"
#include <stdexcept>
#include <sstream>
#include <string.h>

using namespace desres::msys;

ParamTablePtr ParamTable::create() {
    return ParamTablePtr(new ParamTable);
}

ParamTable::ParamTable()
: _nrows(0)
{}

ParamTable::~ParamTable() {
    for (Id i=0; i<_props.size(); i++) {
        if (_props[i].type==StringType) {
            ValueList& vals=_props[i].vals;
            Id j,n = vals.size();
            for (j=0; j<n; j++) free(vals[j].s);
        }
    }
}

Id ParamTable::propIndex(const String& name) const {
    Id i,n = _props.size();
    for (i=0; i<n; i++) {
        if (_props[i].name==name) return i;
    }
    return BadId;
}

void ParamTable::Property::extend() {
    Value v;
    memset(&v, 0, sizeof(v));
    vals.push_back(v);
}

Id ParamTable::addProp( const String& name, ValueType type) {
    Id index = propIndex(name);
    if (!bad(index)) {
        if (propType(index)!=type) {
            throw std::runtime_error("attempt to change type of prop");
        }
        return index;
    }

    _props.push_back(Property());
    Property& prop = _props.back();
    prop.name = name;
    prop.type = type;
    for (Id i=0; i<_nrows; i++) prop.extend();
    return _props.size()-1;
}

Id ParamTable::addParam() {
    for (Id i=0; i<_props.size(); i++) _props[i].extend();
    return _nrows++;
}

void ParamTable::delProp(Id index) {
    if (bad(index)) return;
    if (index>=_props.size()) {
        std::stringstream ss;
        ss << "delProp: no such index " << index;
        throw std::runtime_error(ss.str());
    }
    _props.erase(_props.begin()+index);
}
