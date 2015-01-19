#include "dtrframe.hxx"
#include "dtrutil.hxx"
#include "endianswap.h"
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <set>

using namespace desres::molfile::dtr;

static const uint32_t magic_frame = 0x4445534d;
static const uint32_t align_size  = 8;
static const uint32_t block_size  = 4096;   // FIXME: make this 1 (8?)
static const uint32_t s_version     = 0x00000100;
static const uint32_t s_irosetta    = 0x12345678;
static const float    s_frosetta    = 1234.5;
static const double   s_drosetta    = 1234.5e6;
static const uint32_t s_lrosetta_lo = 0x89abcdef;
static const uint32_t s_lrosetta_hi = 0x01234567;


// 1:1 correspondence with Key::TYPE_ enum 
static const char* typenames[] = {
    "",
    "int32_t",
    "uint32_t",
    "int64_t",
    "uint64_t",
    "float",
    "double",
    "char",
    "unsigned char"
};

static const unsigned elemsizes[] = {
    0,
    4,
    4,
    8,
    8,
    4,
    8,
    1,
    1
};
static const unsigned ntypenames = sizeof(typenames)/sizeof(typenames[0]);

/*!
 * The byte order associated with this machine.  We use
 * 1234 for little endian, 4321 for big endian, and
 * 3412 for the unlikely PDB endianism.
 */

static inline uint32_t machineEndianism() {
#if __BYTE_ORDER == __LITTLE_ENDIAN
    static const uint32_t byteorder = 1234;
#else
#if __BYTE_ORDER == __BIG_ENDIAN
    static const uint32_t byteorder = 4321;
#else
#ifdef PDB_ENDIAN
#if __BYTE_ORDER == __PDB_ENDIAN
    static const uint32_t byteorder = 3412;
#endif
#endif
#endif
#endif
    // If we get a compile error here, then __BYTE_ORDER
    // has an unexpected value.
    return byteorder;
}


/*!
 * See RFC 1146 for Fletcher's Checksum (http://tools.ietf.org/html/rfc1146)
 */

uint32_t 
desres::molfile::dtr::fletcher( const uint16_t *data, unsigned len ) {
    uint32_t sum1 = 0xffff, sum2 = 0xffff;
    int nchunks = len/360;
    len -= 360 * nchunks;
    while (nchunks--) {
        uint32_t i=360;
        sum2 += i*sum1;
        while (i) {
            uint32_t a = *data++;
            uint32_t b = *data++;
            uint32_t c = *data++;
            uint32_t d = *data++;
            uint32_t e = *data++;
            uint32_t f = *data++;
            uint32_t g = *data++;
            uint32_t h = *data++;
            sum1 += a+b+c+d+e+f+g+h;
            sum2 += (i--)*a;
            sum2 += (i--)*b;
            sum2 += (i--)*c;
            sum2 += (i--)*d;
            sum2 += (i--)*e;
            sum2 += (i--)*f;
            sum2 += (i--)*g;
            sum2 += (i--)*h;
        }
        sum1 = (sum1 & 0xffff) + (sum1 >> 16);
        sum2 = (sum2 & 0xffff) + (sum2 >> 16);
    }
    while (len--) {
        sum1 += *data++;
        sum2 += sum1;
    }
    sum1 = (sum1 & 0xffff) + (sum1 >> 16);
    sum2 = (sum2 & 0xffff) + (sum2 >> 16);
    sum1 = (sum1 & 0xffff) + (sum1 >> 16);
    sum2 = (sum2 & 0xffff) + (sum2 >> 16);
    return sum2 << 16 | sum1;
}

namespace {

    struct header_t {
        uint32_t magic;
        uint32_t version;
        uint32_t framesize_lo;  
        uint32_t framesize_hi;  

        uint32_t headersize;    // size of this header
        uint32_t unused0;   
        uint32_t irosetta;      // integer rosetta
        float    frosetta;      // float rosetta

        uint32_t drosetta_lo;
        uint32_t drosetta_hi;

        uint32_t lrosetta_lo;
        uint32_t lrosetta_hi;

        uint32_t endianism;     // endianism of writer
        uint32_t nlabels;       // number of labeled fields
        uint32_t metasize;      // bytes of meta storage
        uint32_t typesize;      // bytes of type storage

        uint32_t labelsize;     // bytes of label storage
        uint32_t scalarsize;    // bytes of scalar storage
        uint64_t fieldsize;     // bytes of field storage

        uint32_t crcsize;       // size of crc field
        uint32_t padding;       // bytes of padding
        uint32_t unused1;
        uint32_t unused2;
    };

    struct meta_t {
        uint32_t typecode;      // index into types
        uint32_t elemsize;
        uint32_t count_lo;
        uint32_t count_hi;
    };

}

/* experimental dtr version2 stuff */
static const uint32_t magic_frame2 = 0xdeadbeef;

namespace {
    /* the version 2 header. */
    struct header2_t {
        uint32_t magic;     /* 0xdeadbeef, native endian */
        uint32_t checksum;  /* threefry checksum, native endian */
    };
    /* following the header:
     *
     * Labels: NULL-delimited, with a zero-length label terminating the list.
     * The first byte in each label encodes the type, with an index into
     * the typenames array above.  More types may be added later.
     *
     * Data: For each field, uint32_t count, followed by data.  No padding.
     *
     * Overhead for a single byte of data with a label of length 1:
     * 8 byte header, 4 bytes for label, 4 bytes for length = 16 bytes.
     */
}

static KeyMap parse_version2(size_t sz, const void* bytes) {
    KeyMap map;
    return map;
}

std::map<std::string, Key> 
desres::molfile::dtr::ParseFrame(size_t sz, const void* data) {
    std::map<std::string,Key> map;

    /* try version 2 header first */
    if (sz<sizeof(header2_t)) {
        DTR_FAILURE("data is too short");
    }
    if (((const header2_t *)data)->magic == magic_frame2) {
        return parse_version2(sz, data);
    }

    // parse header
    header_t header[1];
    if (sz<sizeof(header)) {
        DTR_FAILURE("data is too short");
    }
    memcpy(header, data, sizeof(header));
    convert_ntohl(header);
    if (header->magic != magic_frame) {
        DTR_FAILURE("frame magic number: got " << header->magic
                << " want " << magic_frame);
    }
    uint64_t meta_start = header->headersize;
    uint64_t type_start = meta_start + header->metasize;
    uint64_t label_start = type_start + header->typesize;
    uint64_t scalar_start = label_start + header->labelsize;
    uint64_t field_start = scalar_start + header->scalarsize;
    uint64_t crc_start = field_start + header->fieldsize;
    if (sz < crc_start+4) {
        DTR_FAILURE("frame is too short: need " << crc_start+4 << " got " << sz);
    }
    const char* bytes = (const char *)data;

    // check crc
    uint32_t crc = *reinterpret_cast<const uint32_t*>(bytes+crc_start);
    if (crc!=0) {
        uint32_t frame_crc = fletcher(reinterpret_cast<const uint16_t*>(bytes),
                crc_start/2);
        if (frame_crc != crc) {
            DTR_FAILURE("checksum failure: want " << crc << " got " << frame_crc);
        }
    }
    if (header->nlabels==0) return map;

    // read type names and convert to enum
    std::vector<int> types;
    for (const char* type=bytes+type_start; *type; type+=1+strlen(type)) {
        unsigned i, n = ntypenames;
        for (i=1; i<n; i++) {
            if (!strcmp(type, typenames[i])) {
                types.push_back(i);
                break;
            }
        }
        if (i==n) {
            DTR_FAILURE("Unrecognized typename '" << type << "' in frame");
        }
    }

    // read disk meta
    const meta_t* meta = reinterpret_cast<const meta_t*>(bytes+meta_start);

    // read labels and associated data
    const char* label = bytes+label_start;
    const char* scalars = bytes+scalar_start;
    const char* fields = bytes+field_start;
    for (uint32_t i=0; i<header->nlabels; i++, label+=1+strlen(label)) {
        uint32_t code = ntohl(meta[i].typecode);
        uint32_t elementsize = ntohl(meta[i].elemsize);
        uint32_t count_lo = ntohl(meta[i].count_lo);
        uint32_t count_hi = ntohl(meta[i].count_hi);
        uint64_t count = assemble64(count_lo,count_hi);
        uint64_t nbytes = elementsize*count;
        const char* addr=NULL;
        if (count<=1) {
            addr = scalars;
            scalars += alignInteger(nbytes, align_size);
        } else {
            addr = fields;
            fields += alignInteger(nbytes, align_size);
        }
        const uint32_t this_endian = machineEndianism();
        const uint32_t that_endian = header->endianism;
        bool swap = false;
        if (this_endian!=that_endian) {
            if ((this_endian==1234 && that_endian==4321) ||
                (this_endian==4321 && that_endian==1234)) {
            swap = true;
            } else {
                DTR_FAILURE("Unsupported frame endianism " << that_endian);
            }
        }
        map[label] = Key(addr, count, types.at(code), swap);
    }

    return map;
}


std::string Key::toString() const {
    const char* buf = reinterpret_cast<const char *>(data);
    return std::string(buf, buf+count);
}

bool Key::get(float* buf) const {
    if (!data) return false;
    if (type==TYPE_FLOAT32) {
        if (buf) memcpy(buf, data, count*sizeof(*buf));
    } else if (type==TYPE_FLOAT64) {
        const double* ptr = (const double*)data;
        if (buf) std::copy(ptr, ptr+count, buf);
    } else {
        DTR_FAILURE("Frame data of type " << typenames[type] << " cannot be converted to float");
    }
    if (swap) swap4_unaligned(buf, count);
    return true;
}

bool Key::get(double* buf) const {
    if (!data) return false;
    if (type==TYPE_FLOAT64) {
        if (buf) memcpy(buf, data, count*sizeof(*buf));
    } else if (type==TYPE_FLOAT32) {
        const float* ptr = (const float*)data;
        if (buf) std::copy(ptr, ptr+count, buf);
    } else {
        DTR_FAILURE("Frame data of type " << typenames[type] << " cannot be converted to double");
    }
    if (swap) swap8_unaligned(buf, count);
    return true;
}

bool Key::get(uint32_t* buf) const {
    if (!data) return false;
    if (type==TYPE_UINT32) {
        if (buf) memcpy(buf, data, count*sizeof(*buf));
    } else {
        DTR_FAILURE("Frame data of type " << typenames[type] << " cannot be converted to uint32");
    }
    if (swap) swap4_unaligned(buf, count);
    return true;
}

bool Key::get(int32_t* buf) const {
    if (!data) return false;
    if (type==TYPE_INT32) {
        if (buf) memcpy(buf, data, count*sizeof(*buf));
    } else {
        DTR_FAILURE("Frame data of type " << typenames[type] << " cannot be converted to int32");
    }
    if (swap) swap4_unaligned(buf, count);
    return true;
}

bool Key::get(uint64_t* buf) const {
    if (!data) return false;
    if (type==TYPE_UINT64) {
        if (buf) memcpy(buf, data, count*sizeof(*buf));
    } else {
        DTR_FAILURE("Frame data of type " << typenames[type] << " cannot be converted to uint64");
    }
    if (swap) swap8_unaligned(buf, count);
    return true;
}

bool Key::get(int64_t* buf) const {
    if (!data) return false;
    if (type==TYPE_INT64) {
        if (buf) memcpy(buf, data, count*sizeof(*buf));
    } else {
        DTR_FAILURE("Frame data of type " << typenames[type] << " cannot be converted to int64");
    }
    if (swap) swap8_unaligned(buf, count);
    return true;
}

void Key::set(const char* buf, uint64_t n) {
    data = (const char *)buf;
    count = n;
    type = TYPE_CHAR;
}
void Key::set(const unsigned char* buf, uint64_t n) {
    data = (const unsigned char *)buf;
    count = n;
    type = TYPE_UCHAR;
}
void Key::set(const int32_t* buf, uint64_t n) {
    data = (const char *)buf;
    count = n;
    type = TYPE_INT32;
}
void Key::set(const uint32_t* buf, uint64_t n) {
    data = (const char *)buf;
    count = n;
    type = TYPE_UINT32;
}
void Key::set(const int64_t* buf, uint64_t n) {
    data = (const char *)buf;
    count = n;
    type = TYPE_INT64;
}
void Key::set(const uint64_t* buf, uint64_t n) {
    data = (const char *)buf;
    count = n;
    type = TYPE_UINT64;
}
void Key::set(const float* buf, uint64_t n) {
    data = (const char *)buf;
    count = n;
    type = TYPE_FLOAT32;
}
void Key::set(const double* buf, uint64_t n) {
    data = (const char *)buf;
    count = n;
    type = TYPE_FLOAT64;
}

static unsigned typename_size(KeyMap const& map) {
    std::set<int> types;
    for (KeyMap::const_iterator m=map.begin(), e=map.end(); m!=e; ++m) {
        types.insert(m->second.type);
    }
    unsigned sz=0;
    for (std::set<int>::iterator t=types.begin(), e=types.end(); t!=e; ++t) {
        sz += strlen(typenames[*t])+1;
    }
    return alignInteger(sz+1, align_size);
}

static unsigned label_size(KeyMap const& map) {
    unsigned sz=0;
    for (KeyMap::const_iterator m=map.begin(), e=map.end(); m!=e; ++m) {
        sz += m->first.size()+1;
    }
    return alignInteger(sz+1, align_size);
}

static unsigned scalar_size(KeyMap const& map) {
    unsigned sz=0;
    for (KeyMap::const_iterator m=map.begin(), e=map.end(); m!=e; ++m) {
        Key const& key = m->second;
        if (key.count<=1) {
            sz += alignInteger(elemsizes[key.type] * key.count, align_size);
        }
    }
    return sz;
}

static uint64_t field_size(KeyMap const& map) {
    uint64_t sz=0;
    for (KeyMap::const_iterator m=map.begin(), e=map.end(); m!=e; ++m) {
        Key const& key = m->second;
        if (key.count>1) {
            sz += alignInteger(elemsizes[key.type] * key.count, align_size);
        }
    }
    return sz;
}

size_t desres::molfile::dtr::ConstructFrame(KeyMap const& map, void ** bufptr) {
    if (!bufptr) return 0;
    uint64_t offset_header_block = 0;
    uint64_t size_header_block =
        alignInteger( sizeof(header_t), align_size );

    uint64_t offset_meta_block = offset_header_block + size_header_block;
    uint64_t size_meta_block =
        alignInteger( map.size()*sizeof(meta_t), align_size );

    uint64_t offset_typename_block = offset_meta_block + size_meta_block;
    uint64_t size_typename_block = typename_size(map);

    uint64_t offset_label_block = offset_typename_block + size_typename_block;
    uint64_t size_label_block = label_size(map);

    uint64_t offset_scalar_block = offset_label_block + size_label_block;
    uint64_t size_scalar_block = scalar_size(map);

    uint64_t offset_field_block = offset_scalar_block + size_scalar_block;
    uint64_t size_field_block = field_size(map);

    uint64_t offset_crc_block = offset_field_block + size_field_block;
    uint64_t size_crc_block = sizeof(uint32_t);

    uint64_t offset_padding_block = offset_crc_block + size_crc_block;
    uint64_t size_padding_block =
        alignInteger(offset_padding_block,block_size) - offset_padding_block;

    uint64_t framesize = offset_padding_block + size_padding_block;

    // construct the frame
    *bufptr = realloc(*bufptr, framesize);
    char* base = (char *)*bufptr;
    memset( base, 0, framesize );

    header_t *header = reinterpret_cast<header_t*>(base+offset_header_block);
    meta_t* diskmeta  = reinterpret_cast<meta_t*>(base+offset_meta_block);
    char*       types = reinterpret_cast<char*>(base+offset_typename_block);
    char*       labels    = reinterpret_cast<char*>(base+offset_label_block);
    char*       scalars   = reinterpret_cast<char*>(base+offset_scalar_block);
    char*       fields    = reinterpret_cast<char*>(base+offset_field_block);
    uint32_t*   crc       = reinterpret_cast<uint32_t*>(base+offset_crc_block);

    /*** header ***/
    memset(header,0,sizeof(header_t));
    header->magic = htonl(magic_frame);
    header->version = htonl(s_version);

    header->framesize_lo = htonl(lobytes(framesize));
    header->framesize_hi = htonl(hibytes(framesize));

    header->headersize = htonl(size_header_block);
    header->unused0 = 0;
    uint64_t lrosetta = assemble64(s_lrosetta_lo,s_lrosetta_hi);
    header->irosetta = s_irosetta;
    header->frosetta = s_frosetta;

    header->drosetta_lo = lobytes(s_drosetta);
    header->drosetta_hi = hibytes(s_drosetta);
    header->lrosetta_lo = lobytes(lrosetta);
    header->lrosetta_hi = hibytes(lrosetta);

    header->endianism = htonl(machineEndianism());
    header->nlabels = htonl(map.size());
    header->metasize = htonl(size_meta_block);
    header->typesize = htonl(size_typename_block);

    header->labelsize = htonl(size_label_block);
    header->scalarsize = htonl(size_scalar_block);
    header->fieldsize = assemble64(
            htonl(lobytes(size_field_block)),
            htonl(hibytes(size_field_block)));

    header->crcsize = htonl(size_crc_block);
    header->padding = htonl(size_padding_block);
    header->unused1 = 0;
    header->unused2 = 0;

    std::map<int, unsigned> typemap;
    for (KeyMap::const_iterator m=map.begin(), e=map.end(); m!=e; ++m) {
        if (typemap.find(m->second.type)==typemap.end()) {
            unsigned code=typemap.size();
            typemap[m->second.type]=code;
            const char *name = typenames[m->second.type];
            types=std::copy(name, name+strlen(name), types);
            *types++ = 0;
        }

        diskmeta->typecode = htonl( typemap[m->second.type] );
        diskmeta->elemsize = htonl( elemsizes[m->second.type]);
        diskmeta->count_lo = htonl( lobytes( m->second.count ));
        diskmeta->count_hi = htonl( hibytes( m->second.count ));
        diskmeta++;

        labels=std::copy(m->first.begin(), m->first.end(), labels);
        *labels++ = 0;

        uint64_t nbytes = m->second.count*elemsizes[m->second.type];
        if (m->second.count <= 1) {
            memcpy( scalars, m->second.data, nbytes );
            scalars += alignInteger( nbytes, align_size );
        } else {
            memcpy( fields, m->second.data, nbytes );
            fields += alignInteger( nbytes, align_size );
        }
    }
    *crc = fletcher(reinterpret_cast<uint16_t*>(base),offset_crc_block/2);
    return framesize;
}



