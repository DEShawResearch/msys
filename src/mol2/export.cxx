#include "../mol2.hxx"
#include "elements.hxx"
#include <boost/foreach.hpp>
#include <stdio.h>
#include <math.h>

using namespace desres::msys;

static void export_ct(SystemPtr mol, Id ct, FILE* fd, 
                      Provenance const& provenance, unsigned flags) {

    Id atype = mol->atomPropIndex("sybyl_type");
    Id btype = mol->bondPropIndex("sybyl_type");

    component_t const& cmp = mol->ct(ct);
    IdList atoms = mol->atomsForCt(ct);
    IdList bonds = mol->bondsForCt(ct);
    IdList residues = mol->residuesForCt(ct);
    std::sort(atoms.begin(), atoms.end());
    std::sort(residues.begin(), residues.end());

    /* molecule record */
    fprintf(fd, "@<TRIPOS>MOLECULE\n");
    fprintf(fd, "%s\n", cmp.name().c_str());
    fprintf(fd, " %4lu %4lu %4lu\n", 
            atoms.size(), bonds.size(), residues.size());
    fprintf(fd, "%s\n", residues.size()==1 ? "SMALL" : "BIOPOLYMER");
    fprintf(fd, "USER_CHARGES\n");
    fprintf(fd, "%s\n", "****");    /* status bits */
    fprintf(fd, "%s: %s\n", 
            provenance.version.c_str(),
            provenance.cmdline.c_str());

    fprintf(fd, "@<TRIPOS>ATOM\n");
    for (Id i=0; i<atoms.size(); i++) {
        Id id = atoms[i];
        const atom_t& atm = mol->atom(id);
        const residue_t& res = mol->residue(atm.residue);

        int atomid = 1 + std::lower_bound(
                atoms.begin(), atoms.end(), id)-atoms.begin();
        int resid = 1 + std::lower_bound(
                residues.begin(), residues.end(), atm.residue)-residues.begin();

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
           atomid, aname.c_str(), atm.x, atm.y, atm.z,
           type.c_str(), resid, rname.c_str(), atm.charge);
    }

    /* bond records */
    fprintf(fd, "@<TRIPOS>BOND\n");
    for (Id i=0; i<bonds.size(); i++) {
        Id id = bonds[i];
        bond_t const& bnd = mol->bond(id);
        Id ai = bnd.i;
        Id aj = bnd.j;
        Id si = std::lower_bound(atoms.begin(), atoms.end(), ai)-atoms.begin();
        Id sj = std::lower_bound(atoms.begin(), atoms.end(), aj)-atoms.begin();
        if (si==atoms.size() || sj==atoms.size()) {
            MSYS_FAIL("Ct " << ct << " has bonds which cross ct boundaries.  Cannot export to MOL2.");
        }
        std::string type;
        if (bad(btype)) {
            std::stringstream ss;
            ss << bnd.order;
            type = ss.str();
        } else {
            type = mol->bondPropValue(i,btype).asString();
        }
        fprintf(fd, "%5u %5u %5u %s\n", i+1, si+1, sj+1, type.c_str());
    }

    /* substructure */
    fprintf(fd, "@<TRIPOS>SUBSTRUCTURE\n");
    for (Id i=0; i<residues.size(); i++) {
        Id id = residues[i];
        residue_t const& res = mol->residue(id);
        chain_t const& chn = mol->chain(res.chain);
        Id ri = mol->atomsForResidue(id).at(0);
        Id si = std::lower_bound(atoms.begin(), atoms.end(), ri)-atoms.begin();
        fprintf(fd, "%7u %-4s %7u %-8s 1 %-4s\n",
                i+1,                                /* residue id */
                res.name.c_str(),                   /* residue name */
                si+1,                               /* root atom */
                residues.size()==1 ? "GROUP" : "RESIDUE",
                chn.name.c_str());
    }
}

void desres::msys::ExportMol2( SystemPtr mol, std::string const& path,
                               Provenance const& provenance,
                               unsigned flags) {

    const char* mode = flags & Mol2Export::Append ? "ab" : "wb";
    FILE* fd = fopen(path.c_str(), mode);
    if (!fd) MSYS_FAIL("Could not open '" << "' for writing.");
    boost::shared_ptr<FILE> dtor(fd, fclose);
    for (Id ct=0; ct<mol->ctCount(); ct++) {
        export_ct(mol,ct,fd,provenance,flags);
    }
}

