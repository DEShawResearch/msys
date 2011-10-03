#ifndef msys_dms_hxx
#define msys_dms_hxx

#include "system.hxx"

namespace desres { namespace msys {

    SystemPtr ImportDMS(const std::string& path, bool structure_only=false);

    void ExportDMS(SystemPtr sys, const std::string& path);
}}

#endif
