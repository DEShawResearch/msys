#ifndef desres_msys_smallstring_hxx
#define desres_msys_smallstring_hxx

#include <string>
#include <sstream>
#include <stdexcept>
#include <iostream>
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

        SmallString(std::string const& o) {
            *this = o;
        }

        operator std::string() const {
            return std::string(s, s+n);
        }

        bool operator==(SmallString<N> const& o) const {
            return !strcmp(s, o.s);
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

        size_t size() const {
            return strlen(s);
        }

        bool empty() const {
            return *s == '\0';
        }

        char const& operator[](size_t i) const {
            return s[i];
        }

        template <int M>
        friend std::ostream& operator<<(std::ostream& out,
                                        SmallString<M> const& self);
    };
    template <int M>
    std::ostream& operator<<(std::ostream& out, SmallString<M> const& self) {
        out << self.s;
        return out;
    }
}}

#endif
