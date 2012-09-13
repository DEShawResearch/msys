#ifndef desres_msys_alchemical_hxx
#define desres_msys_alchemical_hxx

#include "system.hxx"

namespace desres { namespace msys {

    SystemPtr MakeAlchemical( SystemPtr A, SystemPtr B,
                              std::vector<IdPair> pairs,
                              bool avoid_alchemical_noop = true);

}}

#endif
