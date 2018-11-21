#include "../sdf.hxx"
#include "../elements.hxx"
#include "../clone.hxx"
#include "../fastjson/print.hxx"
#include <errno.h>
#include <math.h>

using namespace desres::msys;
using desres::msys::fastjson::floatify;

/* Write x in %10.4f format to buffer */
static void format_coord(char* buf, float const x) {
    memset(buf,' ',10); // FIXME: set the whole buffer to whitespace first
#if 0
    int sign = x<0 ? -1 : 1;
    int i, n = x*(sign*10000);
    for (i=9; i>5; --i) {
        div_t qr = div(n,10);
        buf[i] = '0'+qr.rem;
        n = qr.quot;
    }
    buf[5] = '.';
    for (i=4; i>=0; --i) {
        div_t qr = div(n,10);
        buf[i] = '0'+qr.rem;
        n = qr.quot;
        if (n==0) break;
    }
    if (n!=0 || (sign<0 && i==0)) {
        MSYS_FAIL("Input " << x << " too large for sdf coordinate");
    }
    if (sign<0) buf[--i]='-';
#else
    char tmp[32];
    sprintf(tmp, "%10.4f", x);
    memcpy(buf, tmp, strlen(tmp));
#endif
}

/* Write x in %3d format to buffer */
static void format_short(char* buf, short x) {
    memset(buf, ' ', 3);
    div_t qr;
    bool neg=x<0;
    if (neg) x=-x;
    qr = div(x,10);
    buf[2]='0'+qr.rem;
    if (qr.quot!=0) {
        qr = div(qr.quot,10);
        buf[1]='0'+qr.rem;
        if (qr.quot!=0) {
            qr = div(qr.quot,10);
            buf[0]='0'+qr.rem;
        }
    }
    if (neg) {
        if (x<10) {
            buf[1]='-';
        } else {
            buf[0]='-';
        }
    }
}


static inline char* append(char* ptr, const char* buf, size_t sz) {
    memcpy(ptr, buf, sz);
    ptr += sz;
    return ptr;
}

// construct an sdf entry from a System assuming one ct
static std::string format_ct( SystemPtr mol ) {
    Id natoms = mol->maxAtomId();
    Id nbonds = mol->maxBondId();
    if (natoms > 999 || nbonds > 999) {
        MSYS_FAIL("too many atoms (" << natoms <<
                  ") or bonds (" << nbonds <<
                  ") for sdf format");
    }
    char cntsbuf[] = "        0  0  1  0            999 V2000\n";
    char atombuf[] = "     .         .         .     X   0  0  0  0  0  0\n";
    char bondbuf[] = "           0  0  0\n";
    char mchgbuf[] = "M  CHG  1        \n";
    std::map<Id,int> chargemap;
    auto& ct = mol->ct(0);
    std::string name = ct.name();

    /* compute needed size for atoms, bonds, formal charges and M_END */
    unsigned bufsize = 0;
    bufsize += name.size()+3; // 3 newlines
    bufsize += natoms*(sizeof(atombuf)-1);
    bufsize += nbonds*(sizeof(bondbuf)-1);
    bufsize += sizeof(cntsbuf)-1;
    for (unsigned i=0; i<natoms; i++) {
        if (mol->atomFAST(i).formal_charge!=0) bufsize += sizeof(mchgbuf)-1;
    }
    bufsize += 7;
    std::string sdf(bufsize, ' ');
    char* ptr = const_cast<char *>(sdf.data());

    // header
    ptr = append(ptr, name.c_str(), name.size());
    ptr = append(ptr, "\n\n\n", 3);

    // counts
    format_short(cntsbuf  , natoms);
    format_short(cntsbuf+3, nbonds);
    ptr = append(ptr, cntsbuf, sizeof(cntsbuf)-1);

    // atoms
    for (unsigned i=0; i<natoms; i++) {
        auto const& atm = mol->atomFAST(i);
        format_coord(atombuf   , atm.x);
        format_coord(atombuf+10, atm.y);
        format_coord(atombuf+20, atm.z);
        const char* abbr = AbbreviationForElement(atm.atomic_number);
        atombuf[31] = abbr[0];
        atombuf[32] = abbr[1] ? abbr[1] : ' ';
        auto fc = atm.formal_charge;
        if (fc!=0) chargemap[i] = fc;
        format_short(atombuf+39, atm.stereo_parity);
        ptr = append(ptr, atombuf, sizeof(atombuf)-1);
    }

    // bonds
    for (unsigned i=0; i<nbonds; i++) {
        auto const& bnd = mol->bondFAST(i);
        memset(bondbuf,' ',9);
        if (bnd.stereo<0) {
            format_short(bondbuf  , bnd.j+1);
            format_short(bondbuf+3, bnd.i+1);
        } else {
            format_short(bondbuf  , bnd.i+1);
            format_short(bondbuf+3, bnd.j+1);
        }
        format_short(bondbuf+6, bnd.aromatic ? 4 : bnd.order);
        if (bnd.stereo<0) {
            format_short(bondbuf+9, -bnd.stereo);
        } else {
            format_short(bondbuf+9,  bnd.stereo);
        }
        ptr = append(ptr, bondbuf, sizeof(bondbuf)-1);
    }
    
    // charges
    for (auto it : chargemap) {
        format_short(mchgbuf+10, it.first+1);
        format_short(mchgbuf+14, it.second);
        ptr = append(ptr, mchgbuf, sizeof(mchgbuf)-1);
    }
    ptr = append(ptr, "M  END\n", 7);

    // additional data fields
    char floatbuf[32];
    for (auto const& key : ct.keys()) {
        sdf += ">  <";
        sdf += key;
        sdf += ">\n";

        ValueRef v = ct.value(key);
        switch (v.type()) {
        default:
        case StringType: 
            sdf += v.c_str();
            break;
        case IntType: 
            sdf += std::to_string(v.asInt());
            break;
        case FloatType: 
#ifdef __APPLE__
            if (isfinite(v.asFloat())) {
#else
            if (std::isfinite(v.asFloat())) {
#endif
                floatify(v.asFloat(), floatbuf); 
            } else {
                sprintf(floatbuf, "%f", v.asFloat());
            }
            sdf += floatbuf;
            break;
        }
        sdf += "\n\n";
    }

    // end of entry
    sdf += "$$$$\n";
    return sdf;
}

static bool is_single_ct(SystemPtr mol) {
    return mol->ctCount()==1 && 
           mol->atomCount() == mol->maxAtomId() &&
           mol->bondCount() == mol->maxBondId();
}

namespace desres { namespace msys {

    std::string FormatSdf( SystemPtr mol ) {
        if (is_single_ct(mol)) {
            return format_ct(mol);
        }
        std::string sdf;
        for (auto ct : mol->cts()) {
            sdf += format_ct(Clone(mol, mol->atomsForCt(ct)));
        }
        return sdf;
    }

    void ExportSdf( SystemPtr mol, std::string const& path, unsigned flags) {
        const char* mode = flags & SdfExport::Append ? "a" : "w";
        std::shared_ptr<FILE> fp(fopen(path.c_str(), mode), fclose);
        if (!fp) {
            MSYS_FAIL("Error opening " << path << " for writing: " 
                    << strerror(errno));
        }
        auto sdf = FormatSdf(mol);
        if (fwrite(sdf.data(), sdf.size(), 1, fp.get()) != 1) {
            MSYS_FAIL("Error writing " << path << ": " << strerror(errno));
        }
    }
}}

