//
// Incorporated into msys sources for convenience.
//
#pragma once

#if __cplusplus < 201103L
#ifndef __STDC_FORMAT_MACROS // really?  Apparently for some gcc's.
#define __STDC_FORMAT_MACROS
#endif
#include <inttypes.h>
#else
#include <cinttypes>
#endif
#include <vector>
#include <string>
#include <cstddef>
#include <cstring>
#include <cassert>
#include <cstdio>
#include <utility>
#include <arpa/inet.h>  // for htonl

// DOCUMENTATION_BEGIN

// ThreeRoe - A 128-bit non-cryptographic hash function loosely based
// on threefy, and therefore, by inheritance, threefish and skein.
// The public API was initially inspired by Jenkins' Spooky hash, but
// it has evolved significantly.  Basic usage is to construct a
// ThreeRoe object, possibly call Update to add data to it, and call
// Final(), digest() or hexdigest() to obtain the hash.  It's
// header-only code.  Just #include and go.  E.g.,
//
//    #include <ThreeRoe/ThreeRoe.hpp>
//
//    void *data = ...;
//    size_t len = ...;
//    ThreeRoe::result ThreeRoe(data, len).Final();
//
// The result_type, returned by Final() is std::pair<uint64_t, uint64_t>.
//
// The constructor may be called with a single argument, as long
// as that argument has data() and size() methods.  E.g.,
//
//    std::vector<T> vdata;
//    ThreeRoe tr1(vdata);
//    std::string stringdata;
//    ThreeRoe tr2(stringdata);
//    std::array<int, 99> adata;
//    ThreeRoe tr3(adata);
//
// If you prefer the digest as raw bytes, use digest() instead of Final():
//
//    unsigned char hash[16];
//    ThreeRoe(data, len).digest(hash);
//
// If you prefer a std::string of 32 hexadecimal characters,
// then use hexdigest() instead of Final():
//
//    std::string hexdg = ThreeRoe(data, len).hexdigest();
//
// If you don't have all (or even any) of the data at construction
// time, more data can be added with the Update() method, e.g.,
//
//    ThreeRoe hasher;
//    hasher.Update(data, len);
//    hasher.Update(moredata, morelen);
//
// Like the constructor, Update() has overloads for any type
// that supplies data() and size() methods.
//
// Final(), digest() and hexdigest() are all const, so it's
// permissible to call them more than once or to add more data with
// Update between calls.  E.g.,
//
//    auto h2 = hasher.Final();     // same as hash, above
//    hasher.Update(yetmoredata, len); 
//    auto h3 = hasher.Final();   // different from h2
// 
// Calling Init on an existing ThreeRoe object is equivalent to
// constructing a new object.  Init() has the same overloads as the
// constructor.

// Quality:
//  - ThreeRoe passes all of SMHasher's tests.
//  - Jenkins' froggy, running with HLEN=6, LMMM=15, LARRAY=18,
//    THREADS=32, BYTES=400, BITS=5 finds no collisions through
//    count=2^45 (pairs=2^78).

// Performance: 
//   ./spd.cpp reports that ThreeRoe::Update runs at about 0.3 cycles per
//   byte on a 3.06GHz Xeon 5667 (Westmere).  I.e., bulk conversion at
//   10 GiB/s.  SMHasher reports that short strings (0-32 bytes) take
//   from 45 to 80 cycles per hash.

// Portability:
//
//   Some non-standard code and post-C++03 constructs have been creeping in:
//     - inttypes.h and PRIx64.
//     - _USE_BSD and <endian.h>
//     - use of decltype to SFINAE disambiguate some of the fancier
//       overloads.
//
//   Otherwise, ThreeRoe.hpp is portable C++03.  It compiles cleanly
//   with g++ -Wall, but it has not undergone serious portability
//   testing.  The executable, bin/trsum, (analogous to md5sum and
//   cksum) requires C++11 and some home-grown desres libraries
//   (system_error_wrapper) to build, but they're inessential.  It
//   would be easy to strip it back to C++03.
//
//   ThreeRoe uses no "fancy" instructions.  It requires a good
//   compiler with aggressive inlining to get good performance.  It
//   benefits from an architecture with a "rotate bits" instruction,
//   and with ILP for integer ops like xor and add (i.e., Intel).
//
//   It is not x86-specific, but it *is* untested on any other
//   architecture.  In order to avoid surprises, there are #ifdefs
//   that prevent it from building on anything other than x86_64.  It
//   wouldn't be hard to make it endian-agnostic, but I don't have the
//   need, nor do I have the machines to test on.  So that effort is
//   deferred till another day.
//
//   A static constant uint32_t member ThreeRoe::SMHasher_Verifier is
//   set to the value computed by SMHasher's 'VerificationTest'.  (see
//   smhasher/KeysetTest.cpp and ./verifier.cpp for details).  Any
//   changes to the algorithm will result in changes to the output of
//   the verifier, so it can be thought of as an "algorthmic
//   fingerprint".
//
//     John Salmon  Jun 6, 2014
// DOCUMENTATION_END

class ThreeRoe{
    // Everything "unique" to ThreeRoe is in the private methods:
    // mixoneblock() and finish() and the constants.  Everything else
    // is just boilerplate to support the API.

    // The public API, which follows Jenkins' Spooky is below.

    // Rotation constants from threefry.h.  See
    // http://deshawresearch.com/resources_random123.html and also
    // skein.h from
    // http://csrc.nist.gov/groups/ST/hash/sha-3/Round3/documents/Skein_FinalRnd.zip
    enum  {
        /* These are the R_256 constants from the Threefish reference sources
           with names changed to R_64x4... */
        R_64x4_0_0=14, R_64x4_0_1=16,
        R_64x4_1_0=52, R_64x4_1_1=57,
        R_64x4_2_0=23, R_64x4_2_1=40,
        R_64x4_3_0= 5, R_64x4_3_1=37,
        R_64x4_4_0=25, R_64x4_4_1=33,
        R_64x4_5_0=46, R_64x4_5_1=12,
        R_64x4_6_0=58, R_64x4_6_1=22,
        R_64x4_7_0=32, R_64x4_7_1=32
    };

    enum  {
        /*
        // Output from skein_rot_search: (srs64_B64-X1000)
        // Random seed = 1. BlockSize = 128 bits. sampleCnt =  1024. rounds =  8, minHW_or=57
        // Start: Tue Mar  1 10:07:48 2011
        // rMin = 0.136. #0325[*15] [CRC=455A682F. hw_OR=64. cnt=16384. blkSize= 128].format   
        */
        R_64x2_0_0=16,
        R_64x2_1_0=42,
        R_64x2_2_0=12,
        R_64x2_3_0=31,
        R_64x2_4_0=16,
        R_64x2_5_0=32,
        R_64x2_6_0=24,
        R_64x2_7_0=21
        /* 4 rounds: minHW =  4  [  4  4  4  4 ]
        // 5 rounds: minHW =  8  [  8  8  8  8 ]
        // 6 rounds: minHW = 16  [ 16 16 16 16 ]
        // 7 rounds: minHW = 32  [ 32 32 32 32 ]
        // 8 rounds: minHW = 64  [ 64 64 64 64 ]
        // 9 rounds: minHW = 64  [ 64 64 64 64 ]
        //10 rounds: minHW = 64  [ 64 64 64 64 ]
        //11 rounds: minHW = 64  [ 64 64 64 64 ] */
    };

    // RotL_64 - rotate bits left.  This should compile down
    //  to a single rolq instruction on x86_64..
    static uint64_t RotL_64(uint64_t x, unsigned int N)
    {
        return (x << (N & 63)) | (x >> ((64-N) & 63));
    }

    void mixoneblock(const char *ks){
        mixoneblock(ks, state);
    }

    // mixoneblock - inject one block (4x64-bits) of data
    //   into the state, s and do some mixing:
    //
    //   Inject the data into the state with +=.
    //   Do the first round of ThreeFry4x64
    //   reinject, rotated with += again
    //   Mix s[0] into s[3] and s[2] into s[1] with xor
    static void mixoneblock(const char *data, uint64_t s[4]){
#if defined(__x86_64__) || defined(__i386__) // or any other little-endian with permissive alignment requirements...
        const uint64_t* d64 = (const uint64_t *)data;
        uint64_t k0=d64[0], k1=d64[1], k2=d64[2], k3=d64[3];
#else
#error "Only x86 has been tested.  To compile for a different architecture, you'll have to tweak this #if and possibly write some endian-swapping and alignment-adjusting code in mixoneblock"
        // These *untested* lines *might* work on a machine with
        // permissive alignment requirements.  They are meant to
        // initialize k0 to k3 with the same numeric values as the
        // above code does on x86_64.
        //const uint64_t* d64 = (const uint64_t *)data; // alignment???
        //uint64_t k0=le64toh(d64[0]), k1=le64toh(d64[1]), k2=le64toh(d64[2]), k3=le64toh(d64[3]);
#endif
        uint64_t k4=RotL_64(k3,R_64x4_2_0), k5=RotL_64(k2,R_64x4_1_0); 
        uint64_t k6=RotL_64(k1,R_64x4_2_1), k7=RotL_64(k0,R_64x4_1_1); 

        s[0] += k0; s[1] += k1; s[2] += k2; s[3] += k3; 

        s[0] += s[1]; s[1] = RotL_64(s[1],R_64x4_0_0); s[1] ^= s[0]; 
        s[2] += s[3]; s[3] = RotL_64(s[3],R_64x4_0_1); s[3] ^= s[2]; 

        s[0] += k4; s[1] += k5; s[2] += k6; s[3] += k7;
        s[3] ^= s[0];
        s[1] ^= s[2];
    }

    // finish - this is just barely enough to pass the
    // SMHasher tests for collisions of short strings.
    // It's a hybrid of rounds of Threefry2x64 and
    // Threefry4x64, which seems to do *just* enough
    // mixing.  It may be wise to pay a slight penalty
    // in short-string performance to get a little
    // more mixing here...
    static std::pair<uint64_t, uint64_t> 
    finish(uint64_t len, const uint64_t s[4]){
        uint64_t s0 = s[0]+s[2];
        uint64_t s1 = s[1]+s[3];
        uint64_t s2 = s[2];
        uint64_t s3 = s[3]+len;

        s0 += s1; s1 = RotL_64(s1,R_64x4_0_0); s1 ^= s0; 
        s2 += s3; s3 = RotL_64(s3,R_64x4_0_1); s3 ^= s2; 

        s0 += s1; s1 = RotL_64(s1,R_64x2_0_0); s1 ^= s0;
        s0 += s1; s1 = RotL_64(s1,R_64x2_1_0); s1 ^= s0;

        s0 += s3; s1 += s2;

        s0 += s1; s1 = RotL_64(s1,R_64x2_2_0); s1 ^= s0;
        s0 += s1; s1 = RotL_64(s1,R_64x2_3_0); s1 ^= s0;

        s0 += s2; s1 += s3;

        s0 += s1; s1 = RotL_64(s1,R_64x2_0_0); s1 ^= s0;
        s0 += s1; s1 = RotL_64(s1,R_64x2_1_0); s1 ^= s0;

        s0 += s3; s1 += s2;

        s0 += s1; s1 = RotL_64(s1,R_64x2_2_0); s1 ^= s0;
        s0 += s1; s1 = RotL_64(s1,R_64x2_3_0); s1 ^= s0;

        return std::make_pair(s0, s1);
    }

public:
    // The public API was inspired by Jenkins' Spooky... with several
    // "improvements".

    // The methods from here on could easily be abstracted into a
    // common framework for a wide class of hash functions.

    // result_type - a pair of uint64_t.
    typedef std::pair<uint64_t, uint64_t> result_type;

    // ThreeRoe - the constructor takes two optional seed arguments.
    ThreeRoe(uint64_t seed1=0, uint64_t seed2 = 0){
        Init(seed1, seed2);
    }

    // Taking inspiration from python's hashlib, the constructor
    // also takes initial data.
    ThreeRoe(const void *data, size_t n, uint64_t seed1=0, uint64_t seed2=0){
        Init(data, n, seed1, seed2);
    }

    // Init() is only needed to *re*-initialize a ThreeRoe.
    // A newly constructed ThreeRoe has already been Init'ed.
    // (Unlike Jenkins, Init returns *this, and there's an
    // overload that takes initial data. )
    ThreeRoe& Init(uint64_t seed1 = 0, uint64_t seed2 = 0){
        state[0] = seed1;
        state[1] = 0x3243F6A8885A308Dull; // pi<<60
        state[2] = seed2;
        state[3] = 0;
        bytes_deferred = 0;
        len = 0;
        return *this;
    }

    // An overload of Init that also Update's the object with some initial data.
    ThreeRoe& Init(const void *data, size_t n, uint64_t seed1=0, uint64_t seed2=0){
        Init(seed1, seed2);
        return Update(data, n);
    }

    // Update - add N bytes of data to the state.
    //  (unlike Jenkins, we return a reference to this, to
    //  facilitate chaining, e.g., h.Update('he').Update('llo')
    ThreeRoe& Update(const void *data, size_t n){
        const char *cdata = (const char *)data;
        const char *end = cdata+n;
        cdata = catchup(cdata, end);
        if( end-cdata >= ptrdiff_t(sizeof(deferred)) )
            cdata = bulk(cdata, end);
        defer(cdata, end);
        len += n;
        return *this;
    }

    // Too-clever-by-half? - A templated Update method that works for
    // any type, V, with data() and size() members, e.g., std::vector,
    // std::string, std::array.
    template<typename V>
    ThreeRoe& Update(const V& v){
        return Update((void*)v.data(), v.size()*sizeof(*v.data()));
    }    


#if __cplusplus >= 201103
    // Templated constructor and and Init methods that also take
    // an argument of any type, V, with data() and size() members.
    //
    // These methods are limited to C++11 because We have to avoid
    // ambiguity of, e.g., ThreeRoe(1234) - does it match
    // ThreeRoe(seed1=0, seed2=0) or does it match the template?  We
    // resolve the ambiguity by insisting that V have a .data() member
    // function.
    //
    // N.B.  While we *could* do this without using C++11
    // capabilities, at the very least we'd have to roll our own
    // type_traits, which I'm strongly disinclined to do.  It's 2015.
    // Use a C++11 compiler.
    template<typename V, typename _dummy=decltype(std::declval<V>().data())>
    ThreeRoe(const V& v, uint64_t seed1=0, uint64_t seed2 = 0){
        Init<V>(v, seed1, seed2);
    }

    template<typename V, typename _dummy=decltype(std::declval<V>().data())>
    ThreeRoe& Init(const V& v, uint64_t seed1=0, uint64_t seed2=0){
        Init(seed1, seed2);
        return Update<V>(v);
    }
#endif

    // Final() returns the hash, (a result_type) of the data Updated
    // so-far.  Final is const, so it can be called more than once
    // without suprise.  Note that it's a silent error to call
    // Jenkins' Spooky::Final more than once, or to call
    // Spooky::Update after Spooky::Final without an intervening Init.
    // Both are ok with ThreeRoe.
    result_type
    Final() const{
        uint64_t s[] = {state[0], state[1], state[2], state[3]};
        // pad any deferred data with zeros and mix it in
        // with mixoneblock.
        if(bytes_deferred){
            char zeros[CHARS_PER_INBLK] = {};
            memcpy_imm(zeros, deferred, bytes_deferred);
            mixoneblock(zeros, s);
        }
            
        return finish(len, s);
    }

    // N.B.  Our first attempt at digest() returned a std::vector.
    // But for short to medium-length messages (up to O(100 bytes)) ,
    // the overhead of constructing and returning a vector (about 100
    // cycles) dominated the cost of the computing the hash (about 50
    // cycles + 1/4 cycle per byte).

    // digest - write the 16 byte digest to the argument p.
    // The first 8 bytes are the bigendian  bytes of Final().first.
    // The last 8 bytes are the bigendian bytes of Final().second.
    // The choice to use 'bigendian' makes digest() consistent
    // with the use of printf("%016x") in hexdigest(), which is,
    // in turn, consistent with earlier versions of trsum.
    void digest(void *p) const{
        result_type pp = Final();
        union {
            unsigned char c[16];
            uint64_t u64[2];
        }u;
        u.u64[0] = _htobe64(pp.first);
        u.u64[1] = _htobe64(pp.second);
        memcpy(p, &u.c[0], 16);
    }

    // hexdigest - return a 32-byte string, consisting of
    // the hex digits (%016x) of Final().first, followed
    // by the hex digits (%016x) of Final().second.
    std::string hexdigest() const{
        result_type p = Final();
        char buf[33];
        sprintf(buf, "%016" PRIx64 "%016" PRIx64, p.first, p.second);
        return std::string(buf, 32);
    }
        
    // N.B.  We used to have two overloads of Hash128, one of which
    // returned a result_type, and the other of which followed the
    // wacky conventions of Spooky::Hash128.  The latter we're
    // definitely better off without.  The former is syntactic sugar.
    // Call ThreeRoe(...).Final(), instead.

    // Hash64 is from the Spooky API.  Unadulterated syntactic sugar, but
    // mostly harmless.  Note that Spooky's API insists that you
    // provide a seed to Hash64 and Hash32, but in ThreeRoe the seed
    // is optional.
    static uint64_t Hash64(const void *data, size_t len, uint64_t seed=0){
        return ThreeRoe(data, len, seed).Final().first;
    }

    // In ThreeRoe/0.08 we changed the way SMHasher calls ThreeRoe: it
    // calls the new, endian-independent digest() method, rather than
    // doing endian-dependent type-punning of the values returned by
    // Final().  Since digest() is bigendian, but our hardware is
    // little-endian, this changed the value of the
    // 'SMHasher_Verifier' even though we didn't change any
    // "observable output" from ThreeRoe::Final.
    //
    // FWIW, the old verifier value was 0x2527b9f0;
    static const uint32_t SMHasher_Verifier = 0x6CE2839D;

    // _htobe64 - htobe64 is in libc in some BSDs and glibc since 2.9,
    // but it's not standard and there are confusing variants
    // involving header files that may or may not exist.  So we just
    // write our own, using __builtin_bswap64 with gcc and htonl,
    // which is standardized by posix, without.
    //
    // _htobe64 is public so we can call it in verifier.cpp.
    static uint64_t _htobe64(uint64_t x){
#if __GNUC__
    #if __BYTE_ORDER__ == __ORDER_BIG_ENDIAN__
        return x;
    #else
        return __builtin_bswap64(x);
    #endif
#else
        const uint64_t fourtytwo = 42;
        // Is it too much to expect that the compiler will
        // figure this out at compile time?
        if( *(unsigned char*)(&fourtytwo) == fourtytwo ){
            const uint64_t hi = htonl(x);
            const uint64_t lo = htonl(x>>32);
            return (hi<<32) | lo;
        }else{
            return x;
        }
#endif
    }

    // ThreeRoe::Version was tied to the underlying algorithm, not the
    // implementation or the API.  So whenever some observable output
    // changed, we should have also bumped ThreeRoe::Version to an
    // *unused* value.  Unfortunately, we've managed the Version so
    // poorly that it's worse than if we had no versioning at all.
    //
    // So with ThreeRoe/0.09, ThreeRoe::Version is gone.  The
    // algorithm that has been unchanged from 0.04 onward is the only
    // algorithm we will ever call 'ThreeRoe'. Forever.  Promise.  New
    // versions may offer a different API, but if they're called
    // 'ThreeRoe', then they'll compute the same thing.  If we have a
    // great new idea, we'll give it a whole new name, e.g.,
    // FourRoe(?), rather than risk further confusion over exactly
    // what ThreeRoe really means.
    //
    // There's one problem: ThreeRoe/0.03 computes something
    // different.  And ThreeRoe/0.03 is compiled into production codes
    // and has been used as a digest/validator for petabytes of stored
    // data.  In Feb 2016 we uninstalled 0.03 from the DESRES garden,
    // and we copied the 0.03 code into the zendo source tree, giving
    // it the name 'ZendoThreeRoe'.  The hope is that it will remain
    // sequestered inside zendo and never leak out.  What could
    // possibly go wrong? :-(.
    //
    // What happened?  
    //
    // ThreeRoe/0.03 was released on 06 Jun 2014 with a perfectly
    // viable algorithm and:
    //
    //   ThreeRoe/0.03 - ThreeRoe::Version=1
    //
    // A few days later, we realized that the algorithm in 0.03 had
    // the disconcerting property that:
    //
    //     ThreeRoe-0.03(<empty string>) -> (0,0)
    //
    // ThreeRoe/0.04 was released on 11 Jun 2014 with a very
    // slightly different algorithm.  In 0.04,
    //
    //    ThreeRoe-0.04(<empty string>) -> (random-looking bits)
    //
    // We weren't sure things had settled down, so we set
    // ThreeRoe::Version=-1 to indicate uncertainty/instability.
    //
    //   ThreeRoe/0.04 - ThreeRoe::Version=-1
    //   ThreeRoe/0.05 - ThreeRoe::Version=-1
    //   ThreeRoe/0.06 - ThreeRoe::Version=-1
    //   ThreeRoe/0.07 - ThreeRoe::Version=-1
    //
    // In Sep 2015, we felt that things had settled down.  So in 0.08
    // we set ThreeRoe::Version back to a positive value.  BUT WE
    // FORGOT THAT ThreeRoe::Version=1 was already taken by 0.03, and
    // that 0.03 was in production use in Zendo.  So
    //
    //   ThreeRoe/0.08 - ThreeRoe::Version=1
    //
    // At which point, ThreeRoe::Version became useless as a way of
    // identifying the underlying algorithm!  
    //
    // In Feb 2016 we realized our mistake, released 0.09 with a
    // promise never to change the algorithm ever again, no
    // ThreeRoe::Version, and this comment describing why.
private:
    static const size_t WORDS_PER_INBLK = 4; // in uint64_t
    static const size_t CHARS_PER_INBLK = WORDS_PER_INBLK * 8;
    uint64_t state[4];
    uint64_t deferred[WORDS_PER_INBLK];
    unsigned int bytes_deferred;
    uint64_t len;

    // catchup - use the bytes in range(b,e) to catch up with any bytes
    //  that had previously been deferred.  If there is no pending
    //  deferral, do nothing.  If there are enough bytes in the range
    //  to fill the deferral block, then copy that many, and then
    //  consume the now-full deferral block.  Otherwise, copy all
    //  the bytes from the range to the deferral block.
    const char *catchup(const char *b, const char *e){
        if( bytes_deferred ){
            if( e-b >= ptrdiff_t(sizeof(deferred)-bytes_deferred) ){
                size_t ncpy = sizeof(deferred)-bytes_deferred;
                memcpy_imm((char *)deferred + bytes_deferred, b, ncpy);
                mixoneblock((const char *)deferred);
                bytes_deferred = 0;
                return b+ncpy;
            }else{
                defer(b, e);
                return e;
            }
        }
        return b;
    }
                   
    // defer - Copy the range(b, e) to the deferral block.
    //   e-b must fit in the space available, equal to
    //   sizeof(deferred)-bytes_deferred.
    void defer(const char *b, const char *e){
        assert( (e-b) < ptrdiff_t(sizeof(deferred)-bytes_deferred) );
        memcpy_imm((char *)deferred + bytes_deferred, b, e-b);
        bytes_deferred += (e-b);
    }

    // bulk - Do a bulk injection of an initial prefix of the data in
    //  range(b, e).  Return a pointer to the first byte that *was
    //  not* injected.  It is guaranteed that the unconsumed data is
    //  smaller than CHARS_PER_INBLK chars.
    char *bulk(const char *b, const char *e){
        const int UNROLL=8;
        while( b <= e-UNROLL*CHARS_PER_INBLK ){
            if(UNROLL>=1) mixoneblock(b+0*CHARS_PER_INBLK);
            if(UNROLL>=2) mixoneblock(b+1*CHARS_PER_INBLK);
            if(UNROLL>=3) mixoneblock(b+2*CHARS_PER_INBLK);
            if(UNROLL>=4) mixoneblock(b+3*CHARS_PER_INBLK);
            if(UNROLL>=5) mixoneblock(b+4*CHARS_PER_INBLK);
            if(UNROLL>=6) mixoneblock(b+5*CHARS_PER_INBLK);
            if(UNROLL>=7) mixoneblock(b+6*CHARS_PER_INBLK);
            if(UNROLL>=8) mixoneblock(b+7*CHARS_PER_INBLK);
            b += UNROLL*CHARS_PER_INBLK;
        }
        while( b <= e-CHARS_PER_INBLK){
            mixoneblock(b);
            b += CHARS_PER_INBLK;
        }
        return (char *)b;
    }

    // memcpy_imm - just like memcpy, but for a limited range of N.
    //   gcc seems to optimize memcpy a little more aggressively  if
    //   it knows N at compile-time.
    static void *memcpy_imm(void *to, const void *from, int N){
        switch(N){
#define CASE(NN) case NN: return memcpy(to, from, NN);
            CASE(31);
            CASE(30);
            CASE(29);
            CASE(28);
            CASE(27);
            CASE(26);
            CASE(25);
            CASE(24);
            CASE(23);
            CASE(22);
            CASE(21);
            CASE(20);
            CASE(19);
            CASE(18);
            CASE(17);
            CASE(16);
            CASE(15);
            CASE(14);
            CASE(13);
            CASE(12);
            CASE(11);
            CASE(10);
            CASE(9);
            CASE(8);
            CASE(7);
            CASE(6);
            CASE(5);
            CASE(4);
            CASE(3);
            CASE(2);
            CASE(1);
            CASE(0);
        default: return memcpy(to, from, N);
        }
    }
};
