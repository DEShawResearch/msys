#include "../mol2.hxx"
#include "elements.hxx"
#include <boost/foreach.hpp>
#include <stdio.h>
#include <math.h>

using namespace desres::msys;

void desres::msys::ExportMol2( SystemPtr mol, std::string const& path,
                               Provenance const& provenance,
                               unsigned flags) {

    const char* mode = flags & Mol2Export::Append ? "ab" : "wb";
    FILE* fd = fopen(path.c_str(), mode);
    if (!fd) MSYS_FAIL("Could not open '" << "' for writing.");
    boost::shared_ptr<FILE> dtor(fd, fclose);

    Id atype = mol->atomPropIndex("sybyl_type");
    Id btype = mol->bondPropIndex("sybyl_type");

    /* molecule record */
    fprintf(fd, "@<TRIPOS>MOLECULE\n");
    fprintf(fd, "%s\n", mol->name.c_str());
    fprintf(fd, " %4d %4d %4d\n", 
            mol->atomCount(), mol->bondCount(), mol->residueCount());
    fprintf(fd, "%s\n", mol->residueCount()==1 ? "SMALL" : "BIOPOLYMER");
    fprintf(fd, "USER_CHARGES\n");
    fprintf(fd, "%s\n", "****");    /* status bits */
    fprintf(fd, "%s: %s\n", 
            provenance.version.c_str(),
            provenance.cmdline.c_str());

    /* atom records */
    IdList idmap(mol->maxAtomId(), 0);
    Id index=0;

    fprintf(fd, "@<TRIPOS>ATOM\n");
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
        std::string type = 
            bad(atype) ? AbbreviationForElement(atm.atomic_number)
                       : mol->atomPropValue(i,atype).asString().c_str();

        /* write the atom line */
        fprintf(fd, 
           "%7d %-4s      %8.4f  %8.4f  %8.4f %-6s %4d  %4s  %8.4f\n", 
           index, aname.c_str(), atm.x, atm.y, atm.z,
           type.c_str(), atm.residue+1, rname.c_str(), atm.charge);
    }

    /* bond records */
    fprintf(fd, "@<TRIPOS>BOND\n");
    for (Id i=0; i<mol->maxBondId(); i++) {
        if (!mol->hasBond(i)) continue;
        bond_t const& bnd = mol->bond(i);
        Id ai = idmap[bnd.i];
        Id aj = idmap[bnd.j];
        std::string type;
        if (bad(btype)) {
            std::stringstream ss;
            ss << bnd.order;
            type = ss.str();
        } else {
            type = mol->bondPropValue(i,btype).asString();
        }
        fprintf(fd, "%5u %5u %5u %s\n", i+1, ai, aj, type.c_str());
    }

    /* substructure */
    fprintf(fd, "@<TRIPOS>SUBSTRUCTURE\n");
    for (Id i=0; i<mol->maxResidueId(); i++) {
        if (!mol->atomCountForResidue(i)) continue;
        residue_t const& res = mol->residue(i);
        fprintf(fd, "%7u %-4s %7u %-8s\n",
                i+1,                                /* residue id */
                res.name.c_str(),                   /* residue name */
                mol->atomsForResidue(i).at(0)+1,    /* root atom */
                mol->residueCount()==1 ? "GROUP" : "RESIDUE");
    }
}

