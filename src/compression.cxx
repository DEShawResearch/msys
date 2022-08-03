#include "compression.hxx"
#include "types.hxx"
#include <string.h>
#include <zlib.h>

#ifdef MSYS_WITH_ZSTD
#include <zstd.h>
#endif

#ifdef MSYS_WITH_LZ4
#include <lz4frame.h>
#endif

namespace {

static const int CL_DEFAULT = INT_MAX;

// Base compressed streambuf implementation
class compressed_ostreambuf : public std::streambuf {
    protected:
        std::streambuf *m_sbuf;
        int INCHUNK;
        int OUTCHUNK;

    public:
        compressed_ostreambuf(std::streambuf *buf) : m_sbuf(buf) {}
        virtual ~compressed_ostreambuf() = default;

    protected:
        // default compression functions
        virtual int overflow(int ch) override {
            // We always have an extra character left in the in buffer, put it there
            *pptr() = traits_type::to_char_type(ch);
            pbump(1);
            // write out the data, but expecting more data
            return flush(1);
        }

        virtual int sync() override {
            flush(0);
            return m_sbuf->pubsync();
        }

        virtual int flush(bool more) = 0;
};

#ifdef MSYS_WITH_ZSTD
// Stream buffer for initializing istream with
// For ZSTD compression
class zstd_istreambuf : public std::streambuf {
    std::streambuf* m_sbuf;
    ZSTD_inBuffer inb;
    ZSTD_outBuffer outb;
    ZSTD_DStream *zds;

public:
    zstd_istreambuf(std::streambuf *buf, int);
    ~zstd_istreambuf();
protected:
    int underflow() override;
};

class zstd_ostreambuf : public compressed_ostreambuf {
    ZSTD_inBuffer inb;
    ZSTD_outBuffer outb;
    ZSTD_CStream *zcs;

public:
    zstd_ostreambuf(std::streambuf *buf, int compression_level);
    ~zstd_ostreambuf();
private:
    int flush(bool more) override;
};

#endif

// For GZIP
class gzip_istreambuf : public std::streambuf {
    std::streambuf* m_sbuf;
    z_stream strm;
    static const int CHUNK = 16384;
    char in[CHUNK], out[CHUNK];
    bool closed;

public:
    gzip_istreambuf(std::streambuf *buf, int);
    ~gzip_istreambuf();
protected:
    int underflow() override;
private:
    void close();
};

class gzip_ostreambuf : public compressed_ostreambuf {
    z_stream strm;
    static const int CHUNK = 16384;
    char in[CHUNK], out[CHUNK];

public:
    gzip_ostreambuf(std::streambuf *buf, int compression_level);
    ~gzip_ostreambuf();
private:
    int flush(bool more) override;
};

#ifdef MSYS_WITH_LZ4
// For LZ4
class lz4_istreambuf : public std::streambuf {
    std::streambuf* m_sbuf;
    std::vector<char> src;
    std::vector<char> dst;
    size_t src_idx;
    size_t src_have;
    LZ4F_dctx *m_ctx;
    static const int CHUNK = 64*1024;

    public: 
        lz4_istreambuf(std::streambuf *buf, int);
        ~lz4_istreambuf();
    protected:
        int underflow() override;
};

class lz4_ostreambuf : public compressed_ostreambuf {
    std::vector<char> src;
    std::vector<char> dst;
    LZ4F_cctx *m_ctx;
    LZ4F_preferences_t m_prefs;

    public:
        lz4_ostreambuf(std::streambuf *buf, int compression_level);
        ~lz4_ostreambuf();
    private:
        int flush(bool more) override;
        int overflow(int ch) override;
};

#endif

// Straight file - just point this stream to the same buffer
// as the input stream
class file_istream : public std::istream {
    public:
        file_istream(std::istream &stream) {
            rdbuf(stream.rdbuf());
        }
};

#ifdef MSYS_WITH_ZSTD

zstd_istreambuf::zstd_istreambuf(std::streambuf *buf, int) : m_sbuf(buf)
{
    zds = ZSTD_createDStream();
    if (!zds) {
        MSYS_FAIL("Could not init ZSTD dstream");
    }
    inb.src = malloc(ZSTD_DStreamInSize());
    inb.size = ZSTD_DStreamInSize();
    inb.pos = inb.size;   // empty
    outb.dst = malloc(ZSTD_DStreamOutSize());
    outb.size = ZSTD_DStreamOutSize();
    outb.pos = 0;
    if (!inb.src || !outb.dst) {
        MSYS_FAIL("Could not allocate buffer for ZSTD decompression");
    }
}

int
zstd_istreambuf::underflow()
{
    const int chunk = ZSTD_DStreamInSize();

    while (true) {
        // finished input buf? read more
        if (inb.pos == inb.size) {
            inb.size = m_sbuf->sgetn((char *)inb.src, chunk);
            if (inb.size == 0) {
                return traits_type::eof();
            }
            inb.pos = 0;
        }

        outb.pos = 0;
        size_t ret = ZSTD_decompressStream(zds, &outb, &inb);
        if (ZSTD_isError(ret)) {
            MSYS_FAIL("Error decompressing ZSTD stream: " << ZSTD_getErrorName(ret));
        }

        // Some data written?
        if (outb.pos > 0) {
            // Set buffer pointers to our buffer
            setg((char *)outb.dst, (char *)outb.dst, (char *)outb.dst + outb.pos);
            return traits_type::to_int_type(*gptr());
        }
    }
}

zstd_istreambuf::~zstd_istreambuf()
{
    ZSTD_freeDStream(zds);
    zds = NULL;
    free(outb.dst);
    free((void *)inb.src);
}
                                                             

zstd_ostreambuf::zstd_ostreambuf(std::streambuf *buf, int compression_level) : compressed_ostreambuf(buf)
{
    if (compression_level == CL_DEFAULT) {
        compression_level = ZSTD_CLEVEL_DEFAULT;
    }
    zcs = ZSTD_createCStream();
    if (!zcs) {
        MSYS_FAIL("Could not init ZSTD cstream");
    }
    ZSTD_initCStream(zcs, compression_level);

    INCHUNK = ZSTD_CStreamInSize();
    OUTCHUNK = ZSTD_CStreamOutSize();
    inb.src = malloc(INCHUNK);
    inb.size = 0;
    inb.pos = 0;
    outb.dst = malloc(OUTCHUNK);
    outb.size = OUTCHUNK;
    outb.pos = 0;
    if (!inb.src || !outb.dst) {
        MSYS_FAIL("Could not allocate buffer for ZSTD compression");
    }

    // Point input stream to our buffer, leave 1 byte at the end
    setp((char *)inb.src, (char *)inb.src + INCHUNK - 1);
}

zstd_ostreambuf::~zstd_ostreambuf()
{
    flush(0);   // flush out any remaining data

    ZSTD_freeCStream(zcs);
    zcs = NULL;
    free((void *)inb.src);
    free(outb.dst);
}

int
zstd_ostreambuf::flush(bool more) {
    bool done = false;
    int last_char = 0;
    inb.pos = 0;
    inb.size = pptr() - pbase();
    if (inb.size) last_char = traits_type::to_int_type(*(pptr()-1));
    while (inb.pos != inb.size || !done) {
        int ret = ZSTD_compressStream2(zcs, &outb, &inb, more ? ZSTD_e_continue : ZSTD_e_end);
        if (ZSTD_isError(ret)) {
            MSYS_FAIL("Error compressing zstd " << ZSTD_getErrorName(ret));
        }
        // write compressed data to file
        if (m_sbuf->sputn((const char *)outb.dst, outb.pos) != (int)outb.pos) {
            last_char = traits_type::eof();
        }
        outb.pos = 0;

        // finish compressing till everything is flushed out
        done = (more || ret == 0);
    }
    // now the input buffer is consumed, update the pointers
    setp((char *)inb.src, (char *)inb.src + INCHUNK - 1);
    return last_char;
}

#endif  // ZSTD


gzip_istreambuf::gzip_istreambuf(std::streambuf *buf, int) : m_sbuf(buf)
{
    memset(&strm, 0, sizeof(strm));
    /* magic initialization to read gzip files.  It's not in the zlib 
     * manual, or the example code.  
     * http://stackoverflow.com/questions/1838699/how-can-i-decompress-a-gzip-stream-with-zlib
     */
    if (Z_OK!=inflateInit2(&strm, 16+MAX_WBITS)) {
        MSYS_FAIL("Failed initializing zlib stream");
    }
    closed = 0;
}

gzip_istreambuf::~gzip_istreambuf()
{
    close();
}

int
gzip_istreambuf::underflow()
{
    while (true) {
        if (closed) {    // don't process any more data after stream ends
            return traits_type::eof();
        }
        if (strm.avail_in == 0) {
            strm.avail_in = m_sbuf->sgetn(in, CHUNK);
            if (strm.avail_in == 0) {
                return traits_type::eof();
            }
            strm.next_in = (unsigned char *)in;
        }
        strm.avail_out = CHUNK;
        strm.next_out = (unsigned char *)out;
        int rc = inflate(&strm, Z_NO_FLUSH);
        if (rc == Z_STREAM_END) {
            close();
        }
        else if (rc) {
            MSYS_FAIL("Reading zlib stream failed with rc " << rc << ": " << strm.msg);
        }
        int have = CHUNK - strm.avail_out;
        if (have) {
            // Set buffer pointers to our buffer
            setg(out, out, out + have);
            return traits_type::to_int_type(*gptr());
        }
    }
}

void gzip_istreambuf::close() {
    inflateEnd(&strm);
    closed = 1;
}


gzip_ostreambuf::gzip_ostreambuf(std::streambuf *buf, int compression_level) : compressed_ostreambuf(buf)
{
    if (compression_level == CL_DEFAULT) {
        compression_level = Z_DEFAULT_COMPRESSION;
    }
    INCHUNK = 16384;
    memset(&strm, 0, sizeof(strm));
    // options here are defaults, except to turn on gzip encoding
    if (Z_OK != deflateInit2(&strm, compression_level, Z_DEFLATED, 16+MAX_WBITS, 8, Z_DEFAULT_STRATEGY)) {
        MSYS_FAIL("Failed to initialize zlib stream");
    }
    // Set up input buffer, leave a char at the end
    setp(in, in + INCHUNK - 1);
}

gzip_ostreambuf::~gzip_ostreambuf()
{
    flush(0);
    deflateEnd(&strm);
}

int
gzip_ostreambuf::flush(bool more)
{
    int ret = Z_STREAM_END;
    int last_char = 0;
    strm.next_in = (unsigned char *)in;
    strm.avail_in = pptr() - pbase();
    if (strm.avail_in) last_char = traits_type::to_int_type(*(pptr()-1));
    if (strm.avail_in || !more) {
        do {
            strm.avail_out = CHUNK;
            strm.next_out = (unsigned char *)out;
            ret = deflate(&strm, more ? Z_NO_FLUSH : Z_FINISH);
            if (ret == Z_STREAM_ERROR) {
                MSYS_FAIL("Error in zlib compression");
            }
            // write compressed data to file
            if (m_sbuf->sputn(out, CHUNK - strm.avail_out) != CHUNK - strm.avail_out) {
                last_char = traits_type::eof();
            }
            // continue if the output buffer was full
        } while (strm.avail_out == 0);
    }
    if (!more && ret != Z_STREAM_END) {
        MSYS_FAIL("Failed to finish zlib compression");
    }
    // input buffer is now fully available
    setp(in, in + CHUNK - 1);
    return last_char;
}


#ifdef MSYS_WITH_LZ4

lz4_istreambuf::lz4_istreambuf(std::streambuf *buf, int) : m_sbuf(buf), src_idx(0), src_have(0)
{
    if (LZ4F_isError(LZ4F_createDecompressionContext(&m_ctx, LZ4F_VERSION))) {
        MSYS_FAIL("Could not create LZ4 decompression context");
    }
    src.resize(CHUNK);
    dst.resize(CHUNK*2);
}

lz4_istreambuf::~lz4_istreambuf()
{
    LZ4F_freeDecompressionContext(m_ctx);
}

int
lz4_istreambuf::underflow()
{
    while (true) {
        if (src_have == 0) {
            src_have = m_sbuf->sgetn(src.data(), src.size());
            if (src_have == 0) {
                return traits_type::eof();
            }
            src_idx = 0;
        }

        size_t dst_size = dst.size();
        size_t src_size = src_have;
        size_t ret = LZ4F_decompress(m_ctx, dst.data(), &dst_size, src.data() + src_idx, &src_size, NULL);
        if (LZ4F_isError(ret)) {
            MSYS_FAIL("LZ4 decompression failed: " << LZ4F_getErrorName(ret));
        }
        src_idx += src_size;
        src_have -= src_size;

        if (dst_size) {
            setg(dst.data(), dst.data(), dst.data() + dst_size);
            return traits_type::to_int_type(*gptr());
        }
    }
}

lz4_ostreambuf::lz4_ostreambuf(std::streambuf *buf, int compression_level) : compressed_ostreambuf(buf)
{
    if (compression_level == CL_DEFAULT) {
        compression_level = 0;
    }
    if (LZ4F_createCompressionContext(&m_ctx, LZ4F_VERSION) != 0) {
        MSYS_FAIL("Could not init LZ4 compression");
    }
    memset(&m_prefs, 0, sizeof(m_prefs));
    m_prefs.compressionLevel = compression_level;
    INCHUNK = 64*1024;
    OUTCHUNK = LZ4F_compressBound(INCHUNK, &m_prefs);
}

lz4_ostreambuf::~lz4_ostreambuf()
{
    flush(0);
    LZ4F_freeCompressionContext(m_ctx);
}

// LZ4 overrides overflow because we need to start up the compression stream differently
int
lz4_ostreambuf::overflow(int ch)
{
    // If this is the first time, start compression
    if (src.size() == 0) {
        src.resize(INCHUNK);
        dst.resize(OUTCHUNK);
        setp(src.data(), src.data() + src.size() - 1);
        
        size_t ret = LZ4F_compressBegin(m_ctx, dst.data(), dst.size(), &m_prefs);
        if (LZ4F_isError(ret)) {
            MSYS_FAIL("Error starting LZ4 compression: " << LZ4F_getErrorName(ret));
        }
        if (m_sbuf->sputn(dst.data(), ret) != (int)ret) {
            return traits_type::eof();
        }
    }
    // Use the common overflow logic from here on
    return compressed_ostreambuf::overflow(ch);
}

int
lz4_ostreambuf::flush(bool more)
{
    int last_char = 0;
    size_t have = pptr() - pbase();
    if (have) {
        last_char = traits_type::to_int_type(*(pptr()-1));
        size_t ret = LZ4F_compressUpdate(m_ctx, dst.data(), dst.size(), src.data(), have, NULL);
        if (LZ4F_isError(ret)) {
            MSYS_FAIL("Error in LZ4 compression: " << LZ4F_getErrorName(ret));
        }
        if (ret) {
            if (m_sbuf->sputn(dst.data(), ret) != (int)ret) {
                last_char = traits_type::eof();
            }
        }
        setp(src.data(), src.data() + src.size());
    }

    // Finish compression
    if (!more) {
        size_t ret = LZ4F_compressEnd(m_ctx, dst.data(), dst.size(), NULL);
        if (LZ4F_isError(ret)) {
            MSYS_FAIL("Error finishing LZ4 compression: " << LZ4F_getErrorName(ret));
        }
        if (ret) {
            m_sbuf->sputn(dst.data(), ret);
        }
        // reset back to a state where we can restart compression
        src.resize(0);
        setp(NULL, NULL);
    }
    return last_char;
}

#endif   // LZ4

//
// Wrappers for std::[i/o]stream
//
template <class B, class S>
class stream_wrapper : public S {
    std::unique_ptr<B> buf;
    public:
        stream_wrapper(S &stream, int cl=CL_DEFAULT) : buf(new B(stream.rdbuf(), cl)) {
            // change associated buffer to our compressed stream buffer
            this->rdbuf(buf.get());
        }
        ~stream_wrapper() = default;
};

template <class B> using istream_wrapper = stream_wrapper<B, std::istream>;
template <class B> using ostream_wrapper = stream_wrapper<B, std::ostream>;

} // anonymous namespace

namespace desres { namespace msys {

//
// Public API below
//

//
// Create decompression stream
//
std::unique_ptr<std::istream> 
maybe_compressed_istream(std::istream &file)
{
    std::istream *stream;
    unsigned char magic[4];
    file.read((char *)magic, 4);
    file.seekg(0);
    bool is_gzipped = (file.gcount() >= 2 && magic[0]==0x1f && magic[1]==0x8b);
    bool is_zstd = (file.gcount() >= 4 && magic[0] == 0x28 && magic[1] == 0xb5 && magic[2] == 0x2f && magic[3] == 0xfd);
    bool is_lz4  = (file.gcount() >= 4 && magic[0] == 0x04 && magic[1] == 0x22 && magic[2] == 0x4d && magic[3] == 0x18);

    if (is_gzipped) {
        stream = new istream_wrapper<gzip_istreambuf>(file);
    }
    else if (is_zstd) {
#ifdef MSYS_WITH_ZSTD
        stream = new istream_wrapper<zstd_istreambuf>(file);
#else
        MSYS_FAIL("ZSTD not supported - compile with MSYS_WITH_ZSTD");
#endif
    }
    else if (is_lz4) {
#ifdef MSYS_WITH_LZ4
        stream = new istream_wrapper<lz4_istreambuf>(file);
#else
        MSYS_FAIL("LZ4 not supported - compile with MSYS_WITH_LZ4");
#endif
    }
    else {
        stream = new file_istream(file);
    }
    return std::unique_ptr<std::istream>(stream);
}

//
// Create compressed stream, determines which compression to use based on provided extension name
//
std::unique_ptr<std::ostream>
compressed_ostream(std::ostream &file, const std::string &ext, int compression_level)
{
    std::ostream *stream;
    if (ext == "zst") {
#ifdef MSYS_WITH_ZSTD
        stream = new ostream_wrapper<zstd_ostreambuf>(file, compression_level);
#else
        MSYS_FAIL("ZSTD not supported - compile with MSYS_WITH_ZSTD");
#endif
    }
    else if (ext == "lz4") {
#ifdef MSYS_WITH_LZ4
        stream = new ostream_wrapper<lz4_ostreambuf>(file, compression_level);
#else
        MSYS_FAIL("LZ4 not supported - compile with MSYS_WITH_LZ4");
#endif
    }
    else if (ext == "gz") {
        stream = new ostream_wrapper<gzip_ostreambuf>(file, compression_level);
    }
    else {
        MSYS_FAIL("Don't know which compression to use for extension " << ext);
    }
    return std::unique_ptr<std::ostream>(stream);
}

//
// With default compression level
//
std::unique_ptr<std::ostream>
compressed_ostream(std::ostream &file, const std::string &ext)
{
    return compressed_ostream(file, ext, CL_DEFAULT);
}

//
// Return path extension for supported compression methods, or empty string if not supported
//
std::string
compression_extension(const std::string &path)
{
    if (path.rfind(".gz") == path.size()-3) {
        return "gz";
    }
#ifdef MSYS_WITH_ZSTD
    if (path.rfind(".zst") == path.size()-4) {
        return "zst";
    }
#endif
#ifdef MSYS_WITH_LZ4
    if (path.rfind(".lz4") == path.size()-4) {
        return "lz4";
    }
#endif
    return std::string();
}

}}
