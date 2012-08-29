#include "../mol2.hxx"
#include "elements.hxx"
#include <iostream>
#include <iomanip>

using namespace desres::msys;

const char* guess_atom_type(SystemPtr mol, Id id) {
    return "ANY";
}

static void format_float(std::ostream& out, double v) {
    std::stringstream ss;
    if (v>=0) ss << " ";
    ss.precision(9);
    ss << v;
    out.width(9);
    out << ss.str() << " ";
}

void desres::msys::ExportMol2( SystemPtr mol, std::ostream& out,
                               Provenance const& provenance) {

    /* molecule record */
    out << "@<TRIPOS>MOLECULE" << std::endl;
    out << mol->name << std::endl;
    out << mol->atomCount() << " " 
        << mol->bondCount() << " "
        << mol->residueCount() << std::endl;
    out << "BIOPOLYMER" << std::endl;
    out << "USER_CHARGES" << std::endl;
    out << "****" << std::endl;   /* status bits */
    out << std::endl;   /* comment */

    /* atom records */
    out << "@<TRIPOS>ATOM" << std::endl;
    IdList idmap(mol->maxAtomId(), 0);
    Id index=0;

    for (Id i=0; i<mol->maxAtomId(); i++) {
        if (!mol->hasAtom(i)) continue;
        idmap[i] = ++index;
        const atom_t& atm = mol->atom(i);
        const residue_t& res = mol->residue(atm.residue);

        /* Use atom name unless it's invalid */
        std::string aname(atm.name);
        if (aname.size()<1 || aname.size()>7 || !isalpha(aname[0])) { 
            std::stringstream ss;
            ss << AbbreviationForElement(atm.atomic_number) << i;
            aname = ss.str();
        }

        /* Use residue name unless it's invalid */
        std::string rname(res.name);
        if (rname.size()<1 || rname.size()>7 || !isalpha(rname[0])) {
            rname = "UNK";
        }

        /* guess an atom type */
        const char* type = guess_atom_type(mol, i);

        /* write the atom line */
        out << std::setw(7) << index << " ";
        out << std::setw(7) << aname << " ";
        format_float(out, atm.x);
        format_float(out, atm.y);
        format_float(out, atm.z);
        out << type << " ";
        out.width(4);
        out << res.resid << " ";
        out.width(7);
        out << rname << " ";
        format_float(out, atm.charge);
        out << std::endl;
    }

    /* bond records */
    out << "@<TRIPOS>BOND" << std::endl;
    for (Id i=0; i<mol->maxBondId(); i++) {
        if (!mol->hasBond(i)) continue;
        bond_t const& bnd = mol->bond(i);
        Id ai = idmap[bnd.i];
        Id aj = idmap[bnd.j];
        out << std::setw(7) << i+1 << " ";
        out << std::setw(7) << ai << " ";
        out << std::setw(7) << aj << " ";
        out << "un" << std::endl;
    }

}

