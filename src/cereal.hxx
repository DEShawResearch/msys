#ifndef msys_cereal_hxx
#define msys_cereal_hxx

#include "system.hxx"
#include <iostream>

namespace desres { namespace msys {

    SystemPtr ImportCereal(std::istream& in);
    void ExportCereal(SystemPtr mol, std::ostream& out, Provenance const& prov);

}}

#endif
