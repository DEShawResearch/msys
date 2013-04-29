#ifndef desres_msys_dump_hxx
#define desres_msys_dump_hxx

#include "system.hxx"
#include <stdio.h>

namespace desres { namespace msys {

    struct TextExport {
        enum Flags {
            Default     = 0,
            Reorder     = 1 << 0,
            ParamInfo   = 2 << 0

        };
    };

    /* Write out a text representation of the given system, suitable for
     * diffing.
     *
     * A C-style FILE argument is given here, rather than an ostream,
     * because fprintf is much, much faster than writing to streams. */
    void ExportText(SystemPtr mol, FILE* out, unsigned flags=0);

}}

#endif
