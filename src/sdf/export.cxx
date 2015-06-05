#include "sdf.hxx"
#include "elements.hxx"
#include <fastjson/print.hxx>
#include <boost/format.hpp>
#include <boost/math/special_functions/fpclassify.hpp>
#include <fstream>
#include <math.h>
#include <errno.h>

using namespace desres::msys;
using boost::format;
using desres::msys::fastjson::floatify;

/* Write x in %10.4f format to buffer */
static void format_coord(char* buf, float const x) {
    memset(buf,' ',10); // FIXME: set the whole buffer to whitespace first
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
}

/* Write x in %3d format to buffer */
static void format_short(char* buf, unsigned short const x) {
    div_t qr;
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
}

static inline char* append(char* ptr, const char* buf, size_t sz) {
    memcpy(ptr, buf, sz);
    ptr += sz;
    return ptr;
}

std::string desres::msys::FormatSdf( Molecule const& mol ) {
    if (mol.natoms() > 999 || mol.nbonds() > 999) {
        MSYS_FAIL("too many atoms (" << mol.natoms() <<
                  ") or bonds (" << mol.nbonds() <<
                  ") for sdf format");
    }
    char cntsbuf[] = "        0  0  1  0            999 V2000\n";
    char atombuf[] = "     .         .         .     X   0  0  0  0  0  0\n";
    char bondbuf[] = "           0  0  0\n";

    /* compute needed size */
    unsigned bufsize = 0;
    bufsize += mol.name().size()+3; // 3 newlines
    bufsize += mol.natoms()*(sizeof(atombuf)-1);
    bufsize += mol.nbonds()*(sizeof(bondbuf)-1);
    bufsize += sizeof(cntsbuf)-1;
    for (auto const& v : mol.data()) {
        bufsize += v.first.size() + v.second.size() + 7;
    }
    bufsize += 7+5; // M_END + $$$$.
    std::string sdf(bufsize, ' ');
    char* ptr = const_cast<char *>(sdf.data());

    // header
    ptr = append(ptr, mol.name().c_str(), mol.name().size());
    ptr = append(ptr, "\n\n\n", 3);

    // counts
    format_short(cntsbuf  , mol.natoms());
    format_short(cntsbuf+3, mol.nbonds());
    ptr = append(ptr, cntsbuf, sizeof(cntsbuf)-1);

    // atoms
    for (unsigned i=0, n=mol.natoms(); i<n; i++) {
        Molecule::Atom const& atm = mol.atom(i);
        format_coord(atombuf   , atm.x);
        format_coord(atombuf+10, atm.y);
        format_coord(atombuf+20, atm.z);
        const char* abbr = AbbreviationForElement(atm.atomic_number);
        atombuf[31]=abbr[0];
        if (abbr[1]) atombuf[32]=abbr[1];
        auto fc = atm.formal_charge;
        fc=(fc==0 || fc<-3 || fc>3) ? 0 : 4-fc;
        format_short(atombuf+36, fc);
        ptr = append(ptr, atombuf, sizeof(atombuf)-1);
    }

    // bonds
    for (unsigned i=0, n=mol.nbonds(); i<n; i++) {
        Molecule::Bond const& bnd = mol.bond(i);
        memset(bondbuf,' ',9);
        format_short(bondbuf  , bnd.i+1);
        format_short(bondbuf+3, bnd.j+1);
        format_short(bondbuf+6, bnd.order);
        ptr = append(ptr, bondbuf, sizeof(bondbuf)-1);
    }
    ptr = append(ptr, "M  END\n", 7);
    for (auto const& v : mol.data()) {
        ptr = append(ptr, "> <", 3);
        ptr = append(ptr, v.first.c_str(), v.first.size());
        ptr = append(ptr, ">\n", 2);
        ptr = append(ptr, v.second.c_str(), v.second.size());
        ptr = append(ptr, "\n\n", 2);
    }
    ptr = append(ptr, "$$$$\n", 5);
    return sdf;
}

static void export_ct(SystemPtr mol, Id ct, std::ostream& out) {

    if (mol->atomCountForCt(ct) > 999) {
        MSYS_FAIL("SDF export not support for > 999 atoms, have ct " << ct
                << " with " << mol->atomCountForCt(ct) << " atoms.");
    }
    component_t& cmp = mol->ct(ct);
    IdList atoms = mol->atomsForCt(ct);
    IdList bonds = mol->bondsForCt(ct);
    std::sort(atoms.begin(), atoms.end());

    /* header */
    out << cmp.name() << std::endl;
    out << std::endl;
    out << std::endl;

    /* counts */
    out << format("%3d") % atoms.size()
        << format("%3d") % bonds.size()
        << "  0  0  1  0            999 V2000"
        << std::endl;

    /* atoms */
    for (Id i=0; i<atoms.size(); i++) {
        atom_t const& atm = mol->atom(atoms[i]);
        const char* elem = AbbreviationForElement(atm.atomic_number);
        int fc=atm.formal_charge;
        fc=(fc==0 || fc<-3 || fc>3) ? 0 : 4-fc;
        out << format("%10.4f") % atm.x
            << format("%10.4f") % atm.y
            << format("%10.4f ") % atm.z
            << format("%-3s")   % elem
            << " 0"
            << format("%3d") % fc
            << "  0  0  0  0"
            << std::endl;
    }

    /* bonds */
    for (Id i=0; i<bonds.size(); i++) {
        bond_t const& bnd = mol->bond(bonds[i]);
        Id ai = bnd.i;
        Id aj = bnd.j;
        Id si = std::lower_bound(atoms.begin(), atoms.end(), ai)-atoms.begin();
        Id sj = std::lower_bound(atoms.begin(), atoms.end(), aj)-atoms.begin();
        if (si==atoms.size() || sj==atoms.size()) {
            MSYS_FAIL("Ct " << ct << " has bonds which cross ct boundaries.  Cannot export to SDF.");
        }
        int btype = bnd.order;
        if (bnd.resonant_order==1.5) {
            btype=4;
        }
        out << format("%3i") % (si+1)
            << format("%3i") % (sj+1)
            << format("%3i") % btype
            << "  0  0  0"
            << std::endl;
    }
    /* done with molecule section */
    out << "M  END" << std::endl;

    /* write property block */
    std::vector<String> keys = cmp.keys();
    char floatbuf[32];
    for (unsigned i=0; i<keys.size(); i++) {
        out << ">  <" << keys[i] << ">\n";
        ValueRef v = cmp.value(keys[i]);
        switch (v.type()) {
        default:
        case StringType: 
            out << v.asString(); break;
        case IntType: 
            out << v.asInt(); break;
        case FloatType: 
            if (boost::math::isfinite(v.asFloat())) {
                floatify(v.asFloat(), floatbuf); 
            } else {
                sprintf(floatbuf, "%f", v.asFloat());
            }
            out << floatbuf; break;
        }
        out << "\n\n";
    }
    out << "$$$$" << std::endl;
}

namespace desres { namespace msys {

    void ExportSdf( SystemPtr mol, std::string const& path, unsigned flags) {
        if (path.substr(0,6)=="stdout") {
            ExportSdf(mol, std::cout);
            return;
        } 
        std::ios_base::openmode mode = std::ofstream::out;
        if (flags & SdfExport::Append) {
            mode |= std::ofstream::app;
        }
        std::ofstream out(path.c_str(), mode);
        if (!out) {
            MSYS_FAIL("Error opening " << path << " for writing: "
                    << strerror(errno));
        }
        ExportSdf(mol, out);
    }

    void ExportSdf( SystemPtr mol, std::ostream& out ) {
        for (Id ct=0; ct<mol->ctCount(); ct++) export_ct(mol,ct,out);
    }
}}

