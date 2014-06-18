#include "sdf.hxx"
#include "elements.hxx"
#include <fastjson/print.hxx>
#include <boost/format.hpp>
#include <fstream>
#include <math.h>
#include <errno.h>

using namespace desres::msys;
using boost::format;
using desres::fastjson::floatify;

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
            if (isfinite(v.asFloat())) {
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

