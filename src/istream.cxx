#include "istream.hxx"
#include "types.hxx"
#include <string.h>
#include <zlib.h>

using namespace desres::msys;

namespace {

    class file_istream : public istream {
        std::istream& file;
    public:
        file_istream(std::istream& _file) : file(_file) {}
        void read(char *s, std::streamsize n) { file.read(s,n); }
        std::streamsize gcount() const { return file.gcount(); }
    };

    class gzip_istream : public istream {

        std::istream& file;
        z_stream strm;
        static const int CHUNK = 16384;
        char in[CHUNK], out[CHUNK];
        std::streamsize _gcount;

        void init();
        void close();

    public:
        gzip_istream(std::istream& _file);
        ~gzip_istream() { close(); }
        void read(char *s, std::streamsize n);
        std::streamsize gcount() const { return _gcount; }
    };

}

gzip_istream::gzip_istream(std::istream& _file)
: file(_file), _gcount() {
    memset(&strm, 0, sizeof(strm));
    init();
}

void gzip_istream::init() {
    /* magic initialization to read gzip files.  It's not in the zlib 
     * manual, or the example code.  
     * http://stackoverflow.com/questions/1838699/how-can-i-decompress-a-gzip-stream-with-zlib
     */
    if (Z_OK!=inflateInit2(&strm, 16+MAX_WBITS)) {
        MSYS_FAIL("Failed initializing zlib stream");
    }
}

void gzip_istream::read(char *s, std::streamsize n) {

    _gcount=0;
    while (!strm.avail_out) {
        if (!strm.avail_in) {
            if (file.eof()) return;
            file.read((char *)in,CHUNK);
            strm.avail_in = file.gcount();
            if (strm.avail_in==0) {
                /* interesting.  I guess we're done. */
                _gcount = 0;
                return;
            }
            strm.next_in = (unsigned char *)in;
            if (file.fail() && !file.eof()) {
                throw std::runtime_error("Reading gzipped file failed");
            }
        }
        strm.avail_out = CHUNK;
        strm.next_out = (unsigned char *)out;
        int rc = inflate(&strm, Z_NO_FLUSH);
        if (rc==1) {
            close();
            init();
        } else if (rc) {
            MSYS_FAIL("Reading zlib stream failed with rc " << rc << ": " << strm.msg);
        }
        strm.next_out = (unsigned char *)out;
        strm.avail_out = CHUNK - strm.avail_out;
    }
    std::streamsize have = strm.avail_out;
    std::streamsize nread = std::min(n,have);
    memcpy(s, strm.next_out, nread);
    strm.avail_out -= nread;
    strm.next_out += nread;

    _gcount = nread;
}

void gzip_istream::close() {
    inflateEnd(&strm);
}

istream* istream::wrap(std::istream& file) {
    /* check for gzip magic number */
    bool is_gzipped = file.get()==0x1f && file.get()==0x8b;
    file.seekg(0);
    istream* in = NULL;
    if (is_gzipped) {
        in = new gzip_istream(file);
    } else {
        in = new file_istream(file);
    }
    return in;
}

