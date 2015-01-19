#ifndef desres_msys_dtr_util_hxx
#define desres_msys_dtr_util_hxx

#include <netinet/in.h>     // for ntohl
#include <stdint.h>
#include <string.h>
#include <sstream>
#include <stdexcept>

#ifdef _MSC_VER
#define DTR_LOC __FILE__ << ":" << __LINE__ << "\n" << __FUNCSIG__
#else
#define DTR_LOC __FILE__ << ":" << __LINE__ << "\n" << __PRETTY_FUNCTION__
#endif

#define DTR_FAILURE(args) do { \
    std::stringstream ss; \
    ss << args << "\nlocation: " << DTR_LOC; \
    throw std::runtime_error(ss.str()); \
} while (0)


namespace desres { namespace molfile { namespace dtr {

    /*!
     * Extracts the low 32 bits of a 64 bit integer by masking.
     */
    inline uint32_t lobytes(const uint64_t& x) {
      uint32_t mask = 0xffffffff;
      return x & mask;
    }

    /*!
     * Extract the high 32 bits of a 64 bit integer by shifting.
     */
    inline uint32_t hibytes(const uint64_t& x) {
      return x >> 32;
    }

    /*!
     * Extract the low 32 bits of a 64 bit float as an integer.
     */
    inline uint32_t lobytes(const double& x) {
      union {
        uint64_t ival;
        double   dval;
      } u;
      u.dval = x;
      return lobytes(u.ival);
    }

    /*!
     * Extract the high 32 bits of a 64 bit float as an integer.
     */
    inline uint32_t hibytes(const double& x) {
      union {
        uint64_t ival;
        double   dval;
      } u;
      u.dval = x;
      return hibytes(u.ival);
    }

    inline uint64_t assemble64( uint32_t lo, uint32_t hi) {
        uint64_t hi64 = hi;
        return (hi64 << 32) | lo;
    }

    inline double assembleDouble(uint32_t lo, uint32_t hi) {
        union {
            uint64_t ival;
            double   dval;
        } u;
        u.ival = assemble64(lo,hi);
        return u.dval;
    }

    template <typename T>
    inline void convert_ntohl(T* obj) {
        uint32_t* u = reinterpret_cast<uint32_t*>(obj);
        uint32_t i, n = sizeof(T)/sizeof(n);
        for (i=0; i<n; i++) u[i] = ntohl(u[i]);
    }

    inline uint64_t alignInteger( const uint64_t &x, unsigned border) {
        return x + (border - x%border)%border;
    }

}}}

#endif
