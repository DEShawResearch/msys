#include "vdwmap.hxx"
#include <stdexcept>
#include <sstream>
#include <set>
#include <cstdio>
#include <cmath>    // for HUGE_VAL
#include "ff.hxx"

using namespace desres::msys::mae;
using desres::msys::fastjson::Json;

VdwMap::VdwMap( const Json& ffio_ff ) {
    _rule = ffio_ff.get("ffio_comb_rule").as_string("");
    std::set<std::string> functset;
    const Json& vdwtypes = ffio_ff.get("ffio_vdwtypes");
    if (!vdwtypes) {
        return;
    }
    int i,n = vdwtypes.get("__size__").as_int();
    if (!n) {
        return;
    }
    const Json& functs = vdwtypes.get("ffio_funct");
    const Json& names = vdwtypes.get("ffio_name1").valid() ?
        vdwtypes.get("ffio_name1") : vdwtypes.get("ffio_name");
    for (i=0; i<n; i++) { 
        functset.insert(functs.elem(i).as_string());
        VdwParam param;
        for (int j=1;; j++) {
            char colname[32];
            sprintf(colname, "ffio_c%d", j);
            const Json& col = vdwtypes.get(colname);
            if (!col) break;
            param.push_back(col.elem(i).as_float());
        }
        const char * name = names.elem(i).as_string();
        _params[name] = param;
    }
    if (functset.size()!=1) {
        std::stringstream ss;
        ss << "Expected 1 kind of ffio_funct; got " << functset.size();
        throw std::runtime_error(ss.str());
    }
    _funct = *functset.begin();

    const Json& sites = ffio_ff.get("ffio_atoms").valid() ?
        ffio_ff.get("ffio_atoms") : ffio_ff.get("ffio_sites");
    const Json& types = sites.get("ffio_vdwtype");
    const Json& typesB = sites.get("ffio_vdwtypeB");
    const Json& chargeB = sites.get("ffio_chargeB");
#ifdef DESMOND_USE_SCHRODINGER_MMSHARE
    const Json& chargeC = sites.get("ffio_chargeC");
#endif
    if (types.valid()) {
        int i,n = types.size();
        for (i=0; i<n; i++) {
            _vdwnames.push_back( types.elem(i).as_string());
            if (typesB.valid() && typesB.elem(i).kind()==Json::String) {
                _vdwnamesB.push_back(typesB.elem(i).as_string());
            } else {
                _vdwnamesB.push_back("");
            }
            if (chargeB.valid() && chargeB.elem(i).kind()==Json::Float) {
                _chargeB.push_back(chargeB.elem(i).as_float());
            } else {
                _chargeB.push_back(HUGE_VAL);
            }
#ifdef DESMOND_USE_SCHRODINGER_MMSHARE
            if (chargeC.valid() && chargeC.elem(i).kind()==Json::Float) {
                _chargeC.push_back(chargeC.elem(i).as_float());
            } else {
                _chargeC.push_back(HUGE_VAL);
            }
#endif
        }
    }
    to_lower(_funct);
    to_lower(_rule);

    const Json& combined = ffio_ff.get("ffio_vdwtypes_combined");
    if (combined.valid()) {
        const Json& name1=combined.get("ffio_name1");
        const Json& name2=combined.get("ffio_name2"); 
        if (name1.valid()) {
            int i,n = name1.size();
            for (i=0; i<n; i++) {
                const std::string& n1 = name1.elem(i).as_string();
                const std::string& n2 = name2.elem(i).as_string();
                VdwParam p;
                for (int j=1;; j++) {
                    char colname[32];
                    sprintf(colname, "ffio_c%d", j);
                    const Json& col = combined.get(colname);
                    if (!col) break;
                    p.push_back(col.elem(i).as_float());
                }
                _combined[std::make_pair(n1,n2)]=p;
                _combined[std::make_pair(n2,n1)]=p;
            }
        }
    }
}

const VdwParam& VdwMap::param( const std::string& name ) const {
    
    ParamMap::const_iterator i=_params.find(name);
    if (i!=_params.end()) return i->second;
    std::stringstream ss;
    ss << "Missing vdw parameter with name '" << name << "'";
    throw std::runtime_error(ss.str());
}

const VdwParam& VdwMap::param( const VdwType& name1, 
                               const VdwType& name2 ) const {
    CombinedMap::const_iterator i=_combined.find(std::make_pair(name1,name2));
    if (i==_combined.end()) {
        std::stringstream ss;
        ss << "No combined param for '" << name1 << "', '" << name2 << "'";
        throw std::runtime_error(ss.str());
    }
    return i->second;
}

const VdwType& VdwMap::type( int id ) const {
    int n = _vdwnames.size();
    if (id<1 || id>n) {
        std::stringstream ss;
        ss << "illegal site id " << id;
        throw std::runtime_error(ss.str());
    }
    return _vdwnames[id-1];
}
const VdwType& VdwMap::typeB( int id ) const {
    int n = _vdwnamesB.size();
    if (id<1 || id>n) {
        std::stringstream ss;
        ss << "illegal site id " << id;
        throw std::runtime_error(ss.str());
    }
    return _vdwnamesB[id-1];
}

double VdwMap::chargeB( int id ) const {
    int n = _chargeB.size();
    if (id<1 || id>n) {
        std::stringstream ss;
        ss << "illegal site id " << id;
        throw std::runtime_error(ss.str());
    }
    return _chargeB[id-1];
}

#ifdef DESMOND_USE_SCHRODINGER_MMSHARE
double VdwMap::chargeC( int id ) const {
    int n = _chargeC.size();
    if (id<1 || id>n) {
        std::stringstream ss;
        ss << "illegal site id " << id;
        throw std::runtime_error(ss.str());
    }
    return _chargeC[id-1];
}
#endif