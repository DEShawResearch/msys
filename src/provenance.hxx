#ifndef desres_msys_provenance_hxx
#define desres_msys_provenance_hxx

#include "types.hxx"

namespace desres { namespace msys {

    struct Provenance {
        String timestamp;
        String user;
        String workdir;
        String cmdline;

        static Provenance fromArgs(int argc, char *argv[]);
    };

}}


#endif
