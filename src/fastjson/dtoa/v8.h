#ifndef V8_H_
#define V8_H_

#include <cstdio>
#include <cstring>
//#include <cstdint>
#include <stdint.h>
#include <cassert>
#include <cmath>

#define V8_INFINITY HUGE_VAL
namespace OS {
    static inline double nan_value(){return NAN;}
}


#define ASSERT(x) assert(x)
#define UNREACHABLE() assert(0)
#define UNIMPLEMENTED() assert(0)
#define ceiling ceil


// Helper function used by the CHECK_EQ function when given int
// arguments.  Should not be called directly.
static inline void CheckEqualsHelper(int expected, int value) {
  if (expected != value) {
      assert(0);
  }
}


// Helper function used by the CHECK_EQ function when given int64_t
// arguments.  Should not be called directly.
static inline void CheckEqualsHelper(int64_t expected, int64_t value) {
  if (expected != value) {
      assert(0);
  }
}
// Helper function used by the CHECK function when given string
// arguments.  Should not be called directly.
static inline void CheckEqualsHelper(const char* expected,  const char* value) {
  if((expected == NULL && value != NULL) ||
      (expected != NULL && value == NULL) ||
      (expected != NULL && value != NULL && strcmp(expected, value) != 0)) {
      assert(0);
  }
}

// Helper function used by the CHECK function when given floating
// point arguments.  Should not be called directly.
static inline void CheckEqualsHelper(double expected, double value) {
  // Force values to 64 bit memory to truncate 80 bit precision on IA32.
  volatile double* exp = new double[1];
  *exp = expected;
  volatile double* val = new double[1];
  *val = value;
  if (*exp != *val) {
      assert(0);
  }
  delete[] exp;
  delete[] val;
}


#define CHECK(x) assert(x)
#define CHECK_EQ(x,y) CheckEqualsHelper(x,y)


#define CHECK_GT(x,y) assert((x)>(y))
#define CHECK_GE(x,y) assert((x)>=(y))

// The following macro works on both 32 and 64-bit platforms.
// Usage: instead of writing 0x1234567890123456
//      write V8_2PART_UINT64_C(0x12345678,90123456);
#define V8_2PART_UINT64_C(a, b) (((static_cast<uint64_t>(a) << 32) + 0x##b##u))

#define INLINE(header) inline header  __attribute__((always_inline))

// This is inspired by the static assertion facility in boost.  This
// is pretty magical.  If it causes you trouble on a platform you may
// find a fix in the boost code.
template <bool> class StaticAssertion;
template <> class StaticAssertion<true> { };
// This macro joins two tokens.  If one of the tokens is a macro the
// helper call causes it to be resolved before joining.
#define SEMI_STATIC_JOIN(a, b) SEMI_STATIC_JOIN_HELPER(a, b)
#define SEMI_STATIC_JOIN_HELPER(a, b) a##b
// Causes an error during compilation of the condition is not
// statically known to be true.  It is formulated as a typedef so that
// it can be used wherever a typedef can be used.  Beware that this
// actually causes each use to introduce a new defined type with a
// name depending on the source line.
template <int> class StaticAssertionHelper { };
#define STATIC_CHECK(test)                                                  \
  typedef                                                                   \
    StaticAssertionHelper<sizeof(StaticAssertion<static_cast<bool>(test)>)> \
    SEMI_STATIC_JOIN(__StaticAssertTypedef__, __LINE__)

#define STATIC_ASSERT(test)  STATIC_CHECK(test)

#define DISALLOW_COPY_AND_ASSIGN(TypeName)      \
  TypeName(const TypeName&);                    \
  void operator=(const TypeName&)

// The expression ARRAY_SIZE(a) is a compile-time constant of type
// size_t which represents the number of elements of the given
// array. You should only use ARRAY_SIZE on statically allocated
// arrays.
#define ARRAY_SIZE(a)                                   \
  ((sizeof(a) / sizeof(*(a))) /                         \
  static_cast<size_t>(!(sizeof(a) % sizeof(*(a)))))


#include "utils.h"

#endif
