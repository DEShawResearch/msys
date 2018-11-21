/* @COPYRIGHT@ */

#ifndef viparr_destro_dessert_hpp
#define viparr_destro_dessert_hpp

#include <stdexcept>
#include <cstdlib>
#include <cstring>

// compatibility layer for dessert
namespace desres { namespace msys {

  struct dessert : std::runtime_error {
    explicit dessert(const std::string &s) 
    : std::runtime_error(s.c_str()) {}

    dessert(const std::string& s, const char *file, int line, const char *func)
    : std::runtime_error(s.c_str()) {}

    virtual ~dessert() throw() {}
  };

}}

#ifdef _MSC_VER
#define DESSERT_LOC __FILE__, __LINE__, __FUNCSIG__
#elif defined(_GNUC_)
#define DESSERT_LOC __FILE__, __LINE__, __PRETTY_FUNCTION__
#else
#define DESSERT_LOC __FILE__, __LINE__, "(unknown)"
#endif

#endif
