#ifndef desres_msys_dtr_frame_hxx
#define desres_msys_dtr_frame_hxx

#include <map>
#include <string>
#include <stdint.h>

namespace desres { namespace molfile { namespace dtr {

    struct Key {

        enum {
            TYPE_NONE    = 0
          , TYPE_INT32   = 1
          , TYPE_UINT32  = 2
          , TYPE_INT64   = 3
          , TYPE_UINT64  = 4
          , TYPE_FLOAT32 = 5
          , TYPE_FLOAT64 = 6
          , TYPE_CHAR    = 7
          , TYPE_UCHAR   = 8
        };

        const void* data;   // pointer within original data
        uint64_t    count;  // number of elements
        int         type;   // TYPE_FOO
        bool        swap;   // true if frame endianism is different

        Key() : data(), count(), type(), swap() {}
        Key(const void* d, uint64_t cnt, int t, bool s) 
        : data(d), count(cnt), type(t), swap(s) {}

        std::string toString() const;
        bool get(int32_t* buf) const;
        bool get(uint32_t* buf) const;
        bool get(int64_t* buf) const;
        bool get(uint64_t* buf) const;
        bool get(float* buf) const;
        bool get(double* buf) const;
	bool get(char* buf) const;
	bool get(unsigned char* buf) const;

        void set(const char* buf, uint64_t n);
        void set(const unsigned char* buf, uint64_t n);
        void set(const int32_t* buf, uint64_t n);
        void set(const uint32_t* buf, uint64_t n);
        void set(const int64_t* buf, uint64_t n);
        void set(const uint64_t* buf, uint64_t n);
        void set(const float* buf, uint64_t n);
        void set(const double* buf, uint64_t n);
    };

    typedef std::map<std::string, Key> KeyMap;
    KeyMap ParseFrame(size_t sz, const void* data);

    size_t ConstructFrame(KeyMap const& map, void ** bufptr);

    uint32_t fletcher( const uint16_t *data, unsigned len );

}}}


#endif
