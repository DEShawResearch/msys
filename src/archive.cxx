#include "archive.hxx"
#include <fstream>

#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>

#include <boost/serialization/deque.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/set.hpp>
#include <boost/serialization/shared_ptr.hpp>
#include <boost/serialization/variant.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/weak_ptr.hpp>

#include "snappy.hxx"

#include <msys/version.hxx>

using namespace desres::msys;

void desres::msys::SaveArchive(SystemPtr mol, std::ostream& out) {
    boost::iostreams::filtering_ostream f;
    f.push(snappy_compressor());
    f.push(out);
    boost::archive::binary_oarchive ar(f);
    uint64_t hexversion = MSYS_VERSION_HEX;
    ar << hexversion;
    ar << mol;
}

SystemPtr desres::msys::LoadArchive(std::istream& in) {
    boost::iostreams::filtering_istream f;
    f.push(snappy_decompressor());
    f.push(in);
    boost::archive::binary_iarchive ar(f);
    uint64_t hexversion;
    ar >> hexversion;
    if (hexversion != MSYS_VERSION_HEX) {
        MSYS_FAIL("Archive file's msys version " << std::hex << hexversion 
                  << "differs from current msys version " << MSYS_VERSION_HEX);
    }
    SystemPtr mol = System::create();
    ar >> mol;
    return mol;
}

void desres::msys::ExportArchive(SystemPtr mol, std::string const& path) {
    std::ofstream out(path.c_str());
    SaveArchive(mol, out);
}

SystemPtr desres::msys::ImportArchive(std::string const& path) {
    std::ifstream in(path.c_str());
    return LoadArchive(in);
}
