#ifndef desres_msys_compression_hxx
#define desres_msys_compression_hxx

#include <iostream>
#include <string>
#include <memory>

namespace desres { namespace msys {
    // Convert a possibly compressed stream to an uncompressed stream
    // Supports gzip, LZ4 and ZSTD compression
    // For uncompressed streams, the newly returned stream will share
    // the buffer with the input stream (so uncompressed_stream.rdbuf() == file.rdbuf())
    std::unique_ptr<std::istream> maybe_compressed_istream(std::istream &file);

    // Create compressed stream to given file, pick compression based on extension string
    std::unique_ptr<std::ostream> compressed_ostream(std::ostream &file, const std::string &ext);
    //  .. and with a custom compression level (means different things depending on compression algorithm)
    std::unique_ptr<std::ostream> compressed_ostream(std::ostream &file, const std::string &ext, int compresion_level);

    // Determine if the path extension contains a supported compression method, and return the
    // supported extension (suitable for passing to compressed_ostream above),
    // or empty string if compression is not supported based on this filename
    std::string compression_extension(const std::string &path);
}}


#endif
