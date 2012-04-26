#include "ff.hxx"

using namespace desres::msys::mae;

namespace {
    typedef std::map<std::string, const Ffio *> FfioMap;
    FfioMap& registry() {
        static FfioMap r;
        return r;
    }
}

const Ffio * Ffio::get( const std::string& name ) {
    FfioMap::const_iterator i=registry().find(name);
    if (i!=registry().end()) return i->second;
    return NULL;
}

void Ffio::put( const std::string& name, const Ffio * f ) {
    registry()[name]=f;
}

ParamMap::ParamMap( ParamTablePtr params, 
                    const Json& blk,
                    unsigned ncols,
                    const char** maecols ) 
: _params(params) {
    for (unsigned i=0; i<ncols; i++) {
        _columns.push_back(&blk.get(maecols[i]));
    }
}

ParamMap::ParamMap( ParamTablePtr params,
                    const Json& blk,
                    bool alchemical ) 
: _params(params) {
    unsigned i, ncols = params->propCount();
    if (alchemical) {
        if (ncols % 2) throw std::runtime_error( 
                "expected even number of params in alchemical table");
        ncols /= 2;
    }
    for (i=0; i<ncols; i++) {
        std::string col ="ffio_c"; 
        col += (char)('1' + i);
        _columns.push_back(&blk.get(col.c_str()));
    }
    if (alchemical) for (i=0; i<ncols; i++) {
        std::string col ="ffio_c"; 
        col += (char)('1' + i);
        col += "B";
        _columns.push_back(&blk.get(col.c_str()));
    }
}

Id ParamMap::add(int row) {
    Id n = _columns.size();
    ParamList param(n);
    for (Id i=0; i<n; i++) {
        //std::cout << "row " << row << " col " << i << " " << cols[i]->elem(row) << std::endl;
        param[i] = _columns[i]->elem(row).as_float();
    }
    return add(param);
}

Id ParamMap::add(const ParamList& param) {
    int n = param.size();
    ParamDict::const_iterator iter=_pdict.find(param);
    if (iter!=_pdict.end()) return iter->second;
    Id id = _params->addParam();
    for (int i=0; i<n; i++) {
        _params->value(id,i) = param[i];
    }
    _pdict[param]=id;
    return id;
}
