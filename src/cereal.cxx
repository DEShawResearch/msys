#include "cereal.hxx"
#include "msys/version.hxx"
#include "hash.hxx"
#include <thread>
#include <fstream>
#include <boost/iostreams/device/array.hpp>
#include <boost/iostreams/stream.hpp>
#include "MsysThreeRoe.hpp"
#include "compression.hxx"

// Define this before including cereal to make sure thread safety is on
#define CEREAL_THREAD_SAFE  1

#if __has_include (<cereal/archives/binary.hpp>)
#include <cereal/archives/binary.hpp>
#include <cereal/archives/json.hpp>
#include <cereal/types/memory.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/types/string.hpp>
#include <cereal/types/set.hpp>
#include <cereal/types/map.hpp>
#include <cereal/types/deque.hpp>

// Magic number header for MSYS serailized files
#define CEREAL_MAGIC    { 0x43, 0x45, 0x06, 0x13 }

// Bump this version number when anything in the MSYS serialization changes
#define CEREAL_VERSION  "0.1"

namespace desres { namespace msys {

    SystemPtr ImportCereal(std::string const& path) {
        std::ifstream in(path.c_str());
        if (!in) {
            MSYS_FAIL("Unable to open " << path << " for reading");
        }
    
        return ImportCereal(in);
    }

    SystemPtr ImportCereal(std::istream &_in)
    {
        // Compute hash over serialized data
        ThreeRoe hasher;
        std::unique_ptr<std::istream> in (maybe_compressed_istream(_in));

        static const char cereal_magic[] = CEREAL_MAGIC;
        char magic[4];
        in->read(magic, 4);
        if (!in || in->gcount() != 4 || memcmp(magic, cereal_magic, sizeof(cereal_magic) != 0)) {
            MSYS_FAIL("Could not read Cereal magic number");
        }
        hasher.Update(magic, sizeof(cereal_magic));

        auto mol = System::create();
        Provenance prov;
        // Read in the version first
        {
            cereal::BinaryInputArchive iarchive(*in);
            std::string version;
            iarchive(version);
            if (version != CEREAL_VERSION) {
                MSYS_FAIL("Incompatible serialized version");
            }
            hasher.Update(version);
        }

        std::vector<std::thread> threads(mol->serialize_max());
        // separate data locations for each thread to read into
        std::vector<std::string> datas(mol->serialize_max());
        bool err = 0;
        try {
            for (size_t i = 0; i < mol->serialize_max(); i++) {
                // read size, read data in
                size_t size;
                {
                    cereal::BinaryInputArchive iarchive(*in);
                    iarchive(size);
                    // Note - we don't update the hash here, because we need to update
                    // all the data in order of processing, so omit hashing the size value,
                    // since we're effectively hashing the size by hashing the data itself below
                }
                datas[i].resize(size);
                char *s = (char *)datas[i].data();
                in->read(s, size);
                if ((size_t)in->gcount() != size) {
                    MSYS_FAIL("Failed reading " << size << " bytes from cereal stream, got " << in->gcount());
                }
                threads[i] = std::thread([mol, s, size, i, &err] {
                            try {
                                // This is the only stream library I could find which does not copy the
                                // underlying data and presents it as a stream
                                boost::iostreams::stream<boost::iostreams::array_source> stream(s, size);
                                cereal::BinaryInputArchive iarchive(stream);
                                mol->serialize_x(iarchive, i);
                            } catch (...) {
                                err = 1;
                            }
                            });
            }
        }
        catch (...) {
            // Join any leftover threads, or we'll crash
            for (size_t i = 0; i < threads.size(); i++) {
                if (threads[i].joinable()) threads[i].join();
            }
            throw;
        }
        for (size_t i = 0; i < threads.size(); i++) {
            threads[i].join();
            // Update the hash once we have the data available
            hasher.Update(datas[i]);
        }
        if (err) {
            MSYS_FAIL("Could not decode cereal data");
        }

        // Read in the hash and compare
        {
            cereal::BinaryInputArchive iarchive(*in);
            ThreeRoe::result_type hash;
            iarchive(hash.first, hash.second);
            if (hash != hasher.Final()) {
                MSYS_FAIL("Hash in serialized data does not match");
            }
        }

        return mol;
    }

    void ExportCereal(SystemPtr mol, std::string const& path, Provenance const& prov) {
        std::ofstream fout(path);
        std::unique_ptr<std::ostream> cout;

        if (!fout) {
            MSYS_FAIL("Could not open " << path << " for writing");
        }

        std::string ext = compression_extension(path);
        if (ext.size()) {   // compression
            // specify fast compression, to speed compression and decompression at the cost of file size
            cout = (ext == "zst") ? compressed_ostream(fout, ext, -4) : compressed_ostream(fout, ext);
        }
        ExportCereal(mol, cout ? *cout : fout, prov);
    }


    void ExportCereal(SystemPtr mol, std::ostream &out, Provenance const& prov) {
        // Compute a hash over the serialized data
        ThreeRoe hasher;

        static const char cereal_magic[] = CEREAL_MAGIC;
        out.write(cereal_magic, sizeof(cereal_magic));
        hasher.Update(cereal_magic, sizeof(cereal_magic));

        if (!out) {
            MSYS_FAIL("Could not write magic number to cereal out stream");
        }

        // Write out the version first
        {
            cereal::BinaryOutputArchive oarchive(out);
            std::string version(CEREAL_VERSION);
            oarchive(version);
            hasher.Update(version);
        }

        std::vector<std::stringstream> datas(mol->serialize_max());
        std::vector<std::thread> threads(mol->serialize_max());

        try {
            for (size_t i = 0; i < mol->serialize_max(); i++) {
                std::stringstream &ss = datas[i];
                threads[i] = std::thread([mol, &ss, i] {
                        cereal::BinaryOutputArchive oarchive(ss);
                        mol->serialize_x(oarchive, i);
                        });
            }
            for (size_t i = 0; i < mol->serialize_max(); i++) {
                threads[i].join();
                // Write size, then write data
                {
                    cereal::BinaryOutputArchive oarchive(out);
                    size_t size = datas[i].str().size();
                    oarchive(size);
                    // Note that we omit hashing the size value here - see above
                }
                const std::string &data = datas[i].str();
                out << data;
                hasher.Update(data);
            }
        }
        catch (...) {
            // Must join any leftover threads or we crash
            for (size_t i = 0; i < threads.size(); i++) {
                if (threads[i].joinable()) threads[i].join();
            }
            throw;
        }

        // Finally write out the hash
        {
            cereal::BinaryOutputArchive oarchive(out);
            auto hash = hasher.Final();
            oarchive(hash.first, hash.second);
        }
    }
}}

#else

#warning "NO cereal!"
namespace desres { namespace msys {
    SystemPtr ImportCereal(std::string path) {
        MSYS_FAIL("No cereal support in this build");
    }
    SystemPtr ImportCereal(std::istream& in) {
        MSYS_FAIL("No cereal support in this build");
    }

    void ExportCereal(SystemPtr mol, std::string const& path, Provenance const& prov) {
        MSYS_FAIL("No cereal support in this build");
    }
    void ExportCereal(SystemPtr mol, std::ostream& out, Provenance const& prov) {
        MSYS_FAIL("No cereal support in this build");
    }
}}
#endif





