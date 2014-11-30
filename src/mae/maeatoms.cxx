#include "maeatoms.hxx"
#include <cstdio>

using namespace desres::msys;
using namespace desres::msys::mae;
using desres::msys::fastjson::Json;

MaeAtoms::MaeAtoms( const Json& blk) {

    _name = blk.get("__name__").as_string();
    _nrows = blk.get("__size__").as_int();
    const Json& pseudo = blk.get("ffio_index");
    if (pseudo.valid()) {
        _columns.push_back( &pseudo);
        if (_columns[0]->kind()!=Json::Array) {
            FFIO_ERROR("Bad ffio_index column in " << _name);
        }
    }
    for (int i=0; ; i++) {
        char colname[32];
        sprintf(colname, "ffio_a%c", 'i'+i);
        const Json& col = blk.get(colname);
        if (!col) break;
        if (col.kind()!=Json::Array) {
            FFIO_ERROR("Expected array for column " << colname 
                    << " in block " << _name << "; got " << col.kindstr());
        }
        _columns.push_back(&blk.get(colname));
    }
    _ids.resize(_columns.size());
    _ncols = _ids.size();
}

const IdList& MaeAtoms::ids(int row, int nsites) {
    if (row<0 || row>=_nrows) {
        FFIO_ERROR("Invalid row " << row << " for block " << _name);
    }
    if (nsites<0) nsites=_ncols;
    else if (nsites>_ncols) {
        FFIO_ERROR("Invalid nsites " << nsites << " > " << _ncols 
                << " in " << _name);
    }
    _ids.resize(nsites);
    for (int i=0; i<nsites; i++) {
        try {
            _ids[i] = _columns[i]->elem(row).as_int();
        }
        catch (std::exception& e) {
            const char * colname = _columns[i]->get("__name__").as_string();
            FFIO_ERROR("column " << colname << " row " << row 
                    << " is not an int");
        }
    }
    return _ids;
}

