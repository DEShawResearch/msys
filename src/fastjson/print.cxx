/* @COPYRIGHT@ */

#include "print.hxx"
#include <stdexcept>
#include <fstream>
#include <cstdio>

#include "dtoa/v8.h"
#include "dtoa/utils.h"
#include "dtoa/dtoa.h"

using desres::msys::fastjson::Json;
using desres::msys::fastjson::floatify;

static void indent(std::ostream& out, int depth, std::string const& str) {
    for (; depth; --depth) out << str;
}

static void quotify(std::ostream& out, const std::string& raw) {
    out << '"';
    for (std::string::const_iterator p=raw.begin(); p!=raw.end(); ++p) {
        char c = *p;
        switch (c) {
            case '"':   out << "\\\"" ; break;
            case '\\':  out << "\\\\" ; break;
            default:
                if (iscntrl(c)) {
                    switch (c) {
                        case '\b': out << "\\b"; break;
                        case '\f': out << "\\f"; break;
                        case '\n': out << "\\n"; break;
                        case '\r': out << "\\r"; break;
                        case '\t': out << "\\t"; break;
                        default:
                        {
                            char tmp[8];
                            sprintf(tmp, "\\u%04x", c);
                            out << tmp;
                        }
                    }
                } else {
                    out << c;
                }
        }
    }
    out << '"';
}

void desres::msys::fastjson::floatify(double v, char *p) {
    namespace V8=v8::internal;

    char buf[1+V8::kBase10MaximalLength];
    V8::Vector<char> s(buf,sizeof(buf));
    int sign, length, point;
    V8::DoubleToAscii(v, V8::DTOA_SHORTEST, 0,
                  s, &sign, &length, &point);
    if (sign) *p++ = '-';
    /* length is the number of sigfigs.
     * points is the number of digits printed before the decimal.  If point
     * is zero, we add a leading 0 according to the json spec.  If point
     * is negative, there are zeros after the decimal before we get to
     * our sigfigs.  If point is positive, we consume characters from the
     * buffer.  If point is larger than buffer (length), we add zeros before
     * the decimal.  
     */
    if (point<=0) {
        *p++ = '0';
        *p++ = '.';
        for (int i=0; i<-point; i++) *p++ = '0';
        for (int i=0; i<length; i++) *p++ = s[i];

    } else {
        for (int i=0; i<std::min(length,point); i++) *p++ = s[i];
        for (int i=length; i<point; i++) *p++ = '0';
        *p++ = '.';
        for (int i=point; i<length; i++) *p++ = s[i];
        /* ensure integral values don't end with '.', because the Python
         * json module chokes on them for some reason. */
        if (point>=length) *p++ = '0';
    }
    *p++ = '\0';
}

static void rprint_json( const Json &js, std::ostream& out, int depth, 
                         std::string const& istr, std::string const& nstr ) {

    char buf[18];
    switch (js.kind()) {
        case Json::Invalid:
            throw std::invalid_argument("Got json with invalid type");
            break;
        case Json::Int:
            out << js.as_int();
            break;
        case Json::Bool:
            out << (js.as_bool() ? "true" : "false");
            break;

        case Json::Float:
            floatify( js.as_float(), buf ); 
            out << buf;
            break;

        case Json::String:
            quotify(out, js.as_string());
            break;

        case Json::Null:
            out << "null";
            break;

        case Json::Array:
            if (js.size()==0) { 
                out << "[]";

            } else if (js.size()==1) {
                out << "[";
                rprint_json(js.elem(0), out, depth+1, istr, nstr);
                out << "]";
                
            } else {
                out << "[" << nstr;
                for (int i=0; i<js.size(); i++) {
                    indent(out, depth+1, istr);
                    rprint_json(js.elem(i), out, depth+1, istr, nstr);
                    if (i!=js.size()-1) out << ",";
                    out << nstr;
                }
                indent(out, depth, istr);
                out << "]";
            }
            break;

        case Json::Object:
            out << "{" << nstr;
            for (int i=0; i<js.size(); i++) {
                indent(out, depth+1, istr);
                quotify(out,js.key(i));
                out << " : ";
                rprint_json(js.elem(i), out, depth+1, istr, nstr);
                if (i!=js.size()-1) out << ", ";
                out << nstr;
            }
            indent(out, depth, istr);
            out << "}";
            break;
    }
    if (!out) throw std::runtime_error("writing output failed.");
}

namespace desres { namespace msys { namespace fastjson {

    void print_json( std::ostream& out, const Json& js, 
                     std::string const& indent, std::string const& newline  ) {
        if (!out) throw std::runtime_error("bad output stream");
        rprint_json( js, out, 0, indent, newline );
    }
    void print_json( std::ostream& out, const Json& js ) {
        if (!out) throw std::runtime_error("bad output stream");
        rprint_json( js, out, 0, "    ", "\n" );
    }
    void print_json( const char * path, const Json& js ) {
        std::ofstream out(path);
        print_json( out, js );
        out.close();
    }
}}}

std::ostream& operator<<(std::ostream& o, const desres::msys::fastjson::Json& js) {
    print_json( o, js );
    return o;
}

