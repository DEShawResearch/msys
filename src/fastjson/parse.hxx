#ifndef desres_fastjson_parse_hxx
#define desres_fastjson_parse_hxx

#include "json.hxx"
#include <iostream>

namespace desres { namespace msys { namespace fastjson {

    /* parse input stream */
    void parse_json( std::istream& in, Json& js );
    void parse_json( const char * path, Json& js );

}}}

std::istream& operator>>(std::istream& i, desres::msys::fastjson::Json& a);

#endif
