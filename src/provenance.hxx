#ifndef desres_msys_provenance_hxx
#define desres_msys_provenance_hxx

#include "types.hxx"

namespace desres { namespace msys {

    struct Provenance {
        String version;
        String timestamp;
        String user;
        String workdir;
        String cmdline;
        String executable;

        static Provenance fromArgs(int argc, char *argv[]);

        template <typename Ar> void serialize(Ar& a, unsigned) {
            a & version & timestamp & user & workdir & cmdline & executable;
        }
    };

}}


#endif
