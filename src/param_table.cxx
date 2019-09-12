#include "param_table.hxx"
#include <stdexcept>
#include <sstream>
#include <string.h>
#include <cstdio>

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

void ParamTable::Property::update_index() {
    for (Id i=maxIndexId; i<vals.size(); i++) {
        ValueRef r(type, vals[i]);
        index[r].push_back(i);
    }
    maxIndexId = vals.size();
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
            MSYS_FAIL("could not add prop '" << name << "' with type "
                    << type << "; already exists with type " << propType(index));
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
    _paramrefs.push_back(0);
    return _nrows++;
}

void ParamTable::incref(Id p) {
    if (bad(p)) return;
    if (p>=_paramrefs.size()) {
        MSYS_FAIL("Could not incref param " << p << " in ParamTable of size " << _paramrefs.size());
    }
    ++_paramrefs[p];
}

void ParamTable::decref(Id p) {
    if (bad(p)) return;
    if (p>=_paramrefs.size()) {
        MSYS_FAIL("Invalid param " << p << " in ParamTable of size " << _paramrefs.size());
    }
    if (_paramrefs[p]==0) {
        MSYS_FAIL("param " << p << " already has refcount 0");
    }
    --_paramrefs[p];
}

Id ParamTable::duplicate(Id param) {
    if (bad(param)) return addParam();
    if (!hasParam(param)) {
        std::stringstream ss;
        ss << "ParamTable::duplicate: no such param " << param;
        throw std::runtime_error(ss.str());
    }
    Id dst = addParam();
    for (unsigned i=0; i<_props.size(); i++) {
        _props[i].vals[dst] = _props[i].vals[param];
        if (_props[i].type == StringType) {
            char* s = _props[i].vals[param].s;
            if (s) s=strdup(s);
            _props[i].vals[dst].s = s;
        }

    }
    return dst;
}

int ParamTable::compare(Id L, Id R) {
    if (bad(L) && bad(R)) return 0;
    if (!(hasParam(L) && hasParam(R))) {
        throw std::runtime_error("comparison of invalid param ids");
    }
    Id prop, nprops = propCount();
    for (prop=0; prop<nprops; prop++) {
        int c=value(L, prop).compare(value(R, prop));
        if (c<0) return -1;
        if (c>0) return  1;
    }
    return 0;
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

ValueRef ParamTable::value(Id row, String const& name)  { 
    Id col = propIndex(name);
    if (bad(col)) MSYS_FAIL("No such property '" << name << "'");
    return value(row, col);
}

IdList ParamTable::findInt(Id col, Int const& val) {
    Value v;
    v.i=val;
    ValueRef r(IntType, v);
    return findValue(col,r);
}

IdList ParamTable::findFloat(Id col, Float const& val) {
    Value v;
    v.f=val;
    ValueRef r(FloatType, v);
    return findValue(col,r);
}

IdList ParamTable::findString(Id col, String const& val) {
    Value v;
    v.s=const_cast<char *>(val.c_str());
    ValueRef r(StringType, v);
    return findValue(col,r);
}

IdList ParamTable::findValue(Id col, ValueRef const& v) {
    if (col>=_props.size()) {
        MSYS_FAIL("no such column " << col);
    }
    Property& prop = _props[col];
    prop.update_index();
    Property::Index::iterator p = prop.index.find(v);
    if (p==prop.index.end()) return IdList();
    return p->second;
}
