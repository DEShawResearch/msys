#ifndef desres_msys_clone_hxx
#define desres_msys_clone_hxx

#include "system.hxx"

namespace desres { namespace msys {

    // Bud off a new system from the old.  If ShareParams is true,
    // then ParamTables are never copied; they will be shared between
    // the old and new systems.  By default, copies of the ParamTables
    // are made, but ParamTables shared _within_ the old system will
    // also be shared in the new system.
    struct CloneOption { 
        enum Flags { Default     = 0
                   , ShareParams = 1 << 0
                   , UseIndex    = 1 << 1
        };
    };

    SystemPtr Clone( SystemPtr m, IdList const& atoms,
                     CloneOption::Flags flags = CloneOption::Default );

}}

#endif
