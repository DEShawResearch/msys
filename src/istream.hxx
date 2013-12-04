#ifndef desres_msys_istream_hxx
#define desres_msys_istream_hxx

#include <iostream>

namespace desres { namespace msys {

    /* boost::iostreams::gzip_decompressor lies about its gcount(); it returns
     * the size of the input buffer even on a short read.  That seems to
     * confuse the mae parser.  I've thus rolled my own. 
     *
     * Subclassing std::istream or even std::streambuf isn't at all 
     * straightforward, so I've implemented my own thin interface that
     * provides just what the parser needs.
     */
    class istream {
    public:
        virtual ~istream() {}
        virtual std::streamsize read(char* s, std::streamsize n) = 0;

        static istream* wrap(std::istream& file);
    };

}}

#endif
