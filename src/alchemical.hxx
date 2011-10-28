#ifndef desres_msys_alchemical_hxx
#define desres_msys_alchemical_hxx

#include "system.hxx"

namespace desres { namespace msys {

    SystemPtr CreateAlchemical( SystemPtr A, IdList const& amap,
                                SystemPtr B, IdList const& bmap );

}}

#endif
