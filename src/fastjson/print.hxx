#ifndef desres_fastjson_print_hxx
#define desres_fastjson_print_hxx

#include "json.hxx"
#include <iostream>

namespace desres { namespace msys { namespace fastjson {

    /* write to output stream */
    void print_json( std::ostream& out, const Json& js, 
                     std::string const& indent,
                     std::string const& newline);
    void print_json( std::ostream& out, const Json& js);
    void print_json( const char * path, const Json& js);

    /* format a double into p, which must have at least 18 bytes */
    void floatify(double v, char *p);
}}}

std::ostream& operator<<(std::ostream& o, const desres::msys::fastjson::Json& js);

#endif
