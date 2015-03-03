#pragma once

#include <inttypes.h>
#include <cstddef>
#include <cstring>
#include <cassert>
#include <cstdio>
#include <utility>

// DOCUMENTATION_BEGIN

// ThreeRoe - A 128-bit non-cryptographic hash function loosely based
// on threefy, and therefore, by inheritance, threefish and skein.
// The public API follows Jenkins' Spooky hash with a few aesthetic
// changes.  Basic usage is to construct a ThreeRoe object, call
// Update to add data to it, and call Final to obtain the current
// state of the hash.  It's header-only code.  Just #include and go.
// E.g.,
//
//    #include <ThreeRoe/ThreeRoe.hpp>
//
//    ThreeRoe hasher; // constructor takes two optional uint64_t seeds.
//    ...
//    hasher.Update(data, len);
//    ...
//    hasher.Update(moredata, len);
//    ThreeRoe::result_type hash = hasher.Final();
//
// The result_type is std::pair<uint64_t, uint64_t>.
//
// Final() is const, so it's permissible to call it more than once or
// to add more data with Update after a call to Final.  E.g.,
//
//    auto h2 = hasher.Final();     // same as hash, above
//    hasher.Update(yetmoredata, len); 
//    auto h3 = hasher.Final();   // different from h2
// 
// Calling Init on an existing ThreeRoe object is equivalent to
// constructing a new object.  E.g.,
//    hasher.Init();   // also takes two optional uint64_t seeds.
//
// For convenience, there are static member functions, Hash128 and
// Hash64 that encapsulate the constructor, Update and Final into a
// single call:
//
//    uint64_t h = ThreeRoe::Hash64(data, len);
//    auto hpair = ThreeRoe::Hash128(data, len); // returns result_type, i.e., std::pair of uint64_t.

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
//   Except for extensive use of uint64_t from <inttypes.h>, I think
//   ThreeRoe.hpp is portable C++03.  It compiles cleanly with g++
//   -Wall, but it has not undergone serious portability testing.  The
//   executable, bin/trsum, (analogous to md5sum and cksum) requires
//   C++11 to build.
//
//   It uses no "fancy" instructions.  It requires a good compiler
//   with aggressive inlining to get good performance.  It benefits
//   from an architecture with a "rotate bits" instruction, and with
//   ILP for integer ops like xor and add (i.e., Intel).
//
//   It is not x86-specific, but it *is* endian-sensitve.  I.e., the
//   current version would produce different answers on litle- and
//   big-endian machines.  In order to avoid confusion, there are
//   #ifdefs that prevent it from building on anything other than
//   x86_64.  It wouldn't be hard to make it endian-agnostic, but I
//   don't have the need, nor do I have the machines to test on.  So
//   that effort is deferred till another day.
//
//   A static constant uint32_t member ThreeRoe::SMHasher_Verifier is
//   set to the value computed by SMHasher's 'VerificationTest'.  (see
//   ./verifier.cpp for details).  Any changes to the algorithm will
//   result in changes to this function, so it can be thought of as an
//   "algorthmic fingerprint".  There is also a ThreeRoe::Version,
//   currently 1, which should change if and only if the
//   SMHasher_Verifier changes.  I.e., ThreeRoe::Version is an
//   "algorithmic version", and is not tied to any specific
//   implementation.
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
    //   Do the second round of Threefry4x64
    static void mixoneblock(const char *data, uint64_t s[4]){
#if defined(__x86_64__) || defined(__i386__) // or any other bigendian with permissive alignment requirements...
        const uint64_t* d64 = (const uint64_t *)data;
        uint64_t k0=d64[0], k1=d64[1], k2=d64[2], k3=d64[3];
#else
#error "Only x86 has been tested.  To compile for a different architecture, you'll have to tweak this #if and possibly write some endian-swapping and alignment-ignoring code in mixoneblock"
        //uint64_t k0=??, k1=??, k2=??, k3=??;
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
    // The public API follows Jenkins' Spooky... with some
    // "improvements".

    // The methods from here on could easily be abstracted into a
    // common framework for a wide class of hash functions.

    // result_type - a pair of uint64_t.
    typedef std::pair<uint64_t, uint64_t> result_type;

    // ThreeRoe - the constructor takes two optional seed arguments.
    ThreeRoe(uint64_t seed1=0, uint64_t seed2 = 0){
        Init(seed1, seed2);
    }

    // Init() is only needed to *re*-initialize a ThreeRoe.
    // A newly constructed ThreeRoe has already been Init'ed.
    void Init(uint64_t seed1 = 0, uint64_t seed2 = 0){
        state[0] = seed1;
        state[1] = 0x3243F6A8885A308Dull; // pi<<60
        state[2] = seed2;
        state[3] = 0;
        bytes_deferred = 0;
        len = 0;
    }

    // Update - add N bytes of data to the state.
    void Update(const void *data, size_t n){
        const char *cdata = (const char *)data;
        const char *end = cdata+n;
        cdata = catchup(cdata, end);
        if( end-cdata >= ptrdiff_t(sizeof(deferred)) )
            cdata = bulk(cdata, end);
        defer(cdata, end);
        len += n;
    }

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

    // This is how Spooky defines Final.
    void Final(uint64_t *v1, uint64_t *v2 = 0){
        result_type pair = Final();
        *v1 = pair.first;
        if(v2)
            *v2 = pair.second;
    }

    // Hash128 - returns a std::pair result_type:
    static result_type Hash128(const void *data, size_t len, uint64_t seed = 0){
        ThreeRoe me(seed);
        me.Update(data, len);
        return me.Final();
    }

    // Hash128 - Spooky-style overload, using pointer arguments as
    // inputs (seed) and outputs (hash values).  Yuck!
    static void Hash128(const void  *data, size_t len, uint64_t *io1, uint64_t *io2){
        ThreeRoe me(*io1, *io2);
        me.Update(data, len);
        me.Final(io1, io2);
    }

    // Hash64 and Hash32 complete the Spooky API.  Overkill, but
    // mostly harmless.  Note that Spooky's API insists that you
    // provide a seed to Hash64 and Hash32, but in ThreeRoe the seed
    // is optional.
    static uint64_t Hash64(const void *data, size_t len, uint64_t seed=0){
        return Hash128(data, len, seed).first;
    }

    static uint32_t Hash32(const void *data, size_t len, uint32_t seed=0){
        return Hash64(data, len, seed);
    }

    // This version is tied to the underlying algorithm, not the
    // implementation or the API.  So whenever SMHahser_Verifier, (or
    // any other observable output changes, we should also
    // bump the Version.
    static const int Version = -1; // change to positive when it settles down.
    static const uint32_t SMHasher_Verifier = 0x2527b9f0;

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
