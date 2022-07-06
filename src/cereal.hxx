#ifndef msys_cereal_hxx
#define msys_cereal_hxx

#include "system.hxx"
#include <iostream>

namespace desres { namespace msys {

    SystemPtr ImportCereal(std::istream& in);
    SystemPtr ImportCereal(std::string const& path);
    void ExportCereal(SystemPtr mol, std::string const& path, Provenance const& prov);
    void ExportCereal(SystemPtr mol, std::ostream& out, Provenance const& prov);

}}

#endif
