#include "sdf.hxx"
#include "elements.hxx"
#include <boost/format.hpp>

using boost::format;

namespace desres { namespace msys {

    void ExportSDF( SystemPtr mol, std::ostream& out ) {

        if (mol->atomCount() > 999) {
            MSYS_FAIL("SDF export not support for > 999 atoms, have " << 
                    mol->atomCount());
        }

        /* header */
        out << mol->name << std::endl;
        out << std::endl;
        out << std::endl;

        /* counts */
        out << format("%3d") % mol->atomCount()
            << format("%3d") % mol->bondCount()
            << " 0  0  0  0            999 V2000"
            << std::endl;

        /* atoms */
        IdList idmap(mol->maxAtomId(), BadId);
        Id lid = 0;
        for (Id i=0; i<mol->maxAtomId(); i++) {
            if (!mol->hasAtom(i)) continue;
            idmap[i] = ++lid;
            atom_t const& atm = mol->atom(i);
            const char* elem = AbbreviationForElement(atm.atomic_number);
            out << format("%10.4f") % atm.x
                << format("%10.4f") % atm.y
                << format("%10.4f ") % atm.z
                << format("%-3s")   % elem
                << " 0  0  0  0  0  0"
                << std::endl;
        }

        /* bonds */
        for (Id i=0; i<mol->maxBondId(); i++) {
            if (!mol->hasBond(i)) continue;
            bond_t const& bnd = mol->bond(i);
            out << format("%3i") % idmap.at(bnd.i)
                << format("%3i") % idmap.at(bnd.j)
                << format("%3i") % bnd.order
                << "  0  0  0"
                << std::endl;
        }

        /* write property block? */

        /* done with molecule section */
        out << "M  END" << std::endl;

        out << "$$$$" << std::endl;
    }

}}
