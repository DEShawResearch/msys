#include "regex.hxx"
#include <pcre.h>
#include <stdexcept>
#include <sstream>

using namespace desres::msys::atomsel;

Regex::Regex( const std::string& r, int options ) {
    const char * errmsg;
    int offset;
    /* we only want to find complete matches */
    std::string r2("^(");
    r2 += r;
    r2 += ")$";
    pcre * pc = pcre_compile( 
            r2.c_str(), /* regex */
            options,    /* options */
            &errmsg,    /* error message */
            &offset,    /* location of error */
            NULL );     /* table pointer */
    if (!pc) {
        std::stringstream ss;
        ss << "Error compiling regex: " << errmsg;
        throw std::runtime_error(ss.str());
    }
    regex.reset(pc, pcre_free);
}

bool Regex::match( const std::string& pat ) const {
    int rc = pcre_exec( regex.get(),  /* compiled pattern */
            NULL,         /* study wisdom */
            pat.c_str(),  /* target */
            pat.size(),   /* length of target */
            0,            /* start location */
            0,            /* options */
            NULL,         /* return substrings */
            0 );          /* maxlen of substrings */
    return rc>=0;
}
