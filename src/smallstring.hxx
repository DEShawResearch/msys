#ifndef desres_msys_smallstring_hxx
#define desres_msys_smallstring_hxx

#include <string>
#include <sstream>
#include <stdexcept>
#include <string.h>

namespace desres { namespace msys {

    template <int N>
    class SmallString {
        char s[N+1];
        unsigned char n;

        void _copy(const char* p, size_t sz) {
            if (sz>N) {
                std::stringstream ss;
                ss << "Could not assign string of size " << sz << " ('" 
                    << p << " to SmallString<" << N << ">";
                throw std::runtime_error(ss.str());
            }
            n = sz;
            memcpy(s,p,sz);
            s[sz]=0;
        }

    public:
        SmallString() : n() {
            s[0] = 0;
        }

        operator std::string() const {
            return std::string(s, s+n);
        }

        bool operator==(std::string const& o) const {
            return !strcmp(s, o.data());
        }

        SmallString& operator=(std::string const& p) {
            _copy(p.data(), p.size());
            return *this;
        }

        const char* c_str() const {
            return s;
        }

        char const& operator[](size_t i) const {
            return s[i];
        }


    };
}}

#endif
