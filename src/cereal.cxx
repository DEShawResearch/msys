#include "cereal.hxx"

#if __has_include (<cereal/archives/binary.hpp>)
#include <cereal/archives/binary.hpp>
#include <cereal/types/memory.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/types/string.hpp>
#include <cereal/types/set.hpp>
#include <cereal/types/map.hpp>
#include <cereal/types/deque.hpp>

namespace desres { namespace msys {

    SystemPtr ImportCereal(std::istream& in) {
        auto mol = System::create();
        Provenance prov;
        cereal::BinaryInputArchive iarchive(in);
        iarchive(mol, prov);
        mol->addProvenance(prov);
        return mol;
    }

    void ExportCereal(SystemPtr mol, std::ostream& out, Provenance const& prov) {
        cereal::BinaryOutputArchive oarchive(out);
        oarchive(mol, prov);
    }
}}

#else

#warning "NO cereal!"
namespace desres { namespace msys {
    SystemPtr ImportCereal(std::istream& in) {
        MSYS_FAIL("No cereal support in this build");
    }

    void ExportCereal(SystemPtr mol, std::ostream& out, Provenance const& prov) {
        MSYS_FAIL("No cereal support in this build");
    }
}}
#endif





