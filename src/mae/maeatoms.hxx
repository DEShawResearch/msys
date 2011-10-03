#ifndef desres_msys_mae_maeatoms_hxx
#define desres_msys_mae_maeatoms_hxx

#include "ff.hxx"

namespace desres { namespace msys { namespace mae {

    /* utility class for extracting ids from a forcefield block */
    class MaeAtoms {

        std::vector<const fastjson::Json *> _columns;
        IdList _ids;
        const char * _name;
        int _nrows;
        int _ncols;

    public:
        explicit MaeAtoms(const fastjson::Json& blk);
        const IdList& ids(int row, int nsites=-1);
    };

}}}

#endif
