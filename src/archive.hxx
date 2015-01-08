#ifndef desres_msys_archive_hxx
#define desres_msys_archive_hxx

#include "system.hxx"
#include <iostream>

namespace desres { namespace msys {

    void SaveArchive(SystemPtr mol, std::ostream& out);
    SystemPtr LoadArchive(std::istream& in);

    void ExportArchive(SystemPtr mol, std::string const& path);
    SystemPtr ImportArchive(std::string const& path);
}}

#endif
