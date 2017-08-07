#include "../mol2.hxx"
#include "../elements.hxx"
#include "../clone.hxx"
#include "../annotated_system.hxx"
#include <stdio.h>
#include <math.h>

using namespace desres::msys;

static void export_mol2(SystemPtr omol, AnnotatedSystem& asys,
                        IdList parent_ids,
                        FILE* fd, Provenance const& provenance) {

    if (parent_ids.empty()) {
        parent_ids = omol->atoms();
    }
    SystemPtr mol = Clone(omol, parent_ids);
    std::vector<const char*> atypes(mol->maxAtomId());

    /* molecule record */
    fprintf(fd, "@<TRIPOS>MOLECULE\n");
    fprintf(fd, "%s\n", mol->ct(0).name().data());
    fprintf(fd, " %4u %4u %4u\n", 
            mol->atomCount(), mol->bondCount(), mol->residueCount());
    fprintf(fd, "%s\n", mol->residueCount()==1 ? "SMALL" : "BIOPOLYMER");
    fprintf(fd, "USER_CHARGES\n");
    fprintf(fd, "%s\n", "****");    /* status bits */
    fprintf(fd, "%s: %s\n", 
            provenance.version.data(),
            provenance.cmdline.data());

    /* atoms */
    fprintf(fd, "@<TRIPOS>ATOM\n");
    for (Id local_id=0, n=parent_ids.size(); local_id<n; local_id++) {
        Id parent_id = parent_ids[local_id];
        const atom_t& atm = mol->atomFAST(local_id);
        const residue_t& res = mol->residueFAST(atm.residue);

        /* Use atom name unless it's invalid */
        std::string aname(atm.name);
        if (aname.size()<1 || aname.size()>7 || !isalpha(aname[0])) { 
            aname = AbbreviationForElement(atm.atomic_number)+std::to_string(local_id);
        }

        /* Use residue name unless it's invalid */
        std::string rname(res.name);
        if (rname.size()<1 || rname.size()>7 || !isalpha(rname[0])) {
            rname = "UNK";
        }

        /* guess an atom type using the original molecule */
        bool cyclic = asys.atomRingCount(parent_id);
        const char* type = GuessSybylAtomType(omol, parent_id, cyclic);
        atypes[local_id] = type;

        /* write the atom line */
        const char* bb = atm.type==AtomProBack ? "BACKBONE" : "";
        fprintf(fd, 
           "%4d  %-4s  %8.4f  %8.4f  %8.4f %-5s %4d  %4s%-4d %7.4f %s\n", 
           local_id+1, aname.c_str(), atm.x, atm.y, atm.z,
           type, atm.residue+1, rname.c_str(), res.resid, atm.charge, bb);
    }

    /* bond records */
    fprintf(fd, "@<TRIPOS>BOND\n");
    for (Id local_id = 0, n=mol->bondCount(); local_id<n; local_id++) {
        bond_t const& bnd = mol->bondFAST(local_id);
        Id ai = bnd.i;
        Id aj = bnd.j;
        Id parent_id = omol->findBond(parent_ids[ai], parent_ids[aj]);
        const char* type = GuessSybylBondType(omol, parent_id, atypes[ai], atypes[aj]);
        fprintf(fd, "%5u %5u %5u %s\n", local_id+1, ai+1, aj+1, type);
    }

    /* substructure */
    fprintf(fd, "@<TRIPOS>SUBSTRUCTURE\n");
    for (Id local_id=0, n=mol->residueCount(); local_id<n; local_id++) {
        residue_t const& res = mol->residue(local_id);
        chain_t const& chn = mol->chain(res.chain);
        Id ri = mol->atomsForResidue(local_id).at(0);
        fprintf(fd, "%7u %4s%-4d %4u %-8s 1 %-4s\n",
                local_id+1,                         /* residue id */
                res.name.c_str(), res.resid,        /* residue name */
                ri+1,                               /* root atom */
                n==1 ? "GROUP" : "RESIDUE",
                chn.name.c_str());
    }

    /* additional properties */
    auto keys = mol->ct(0).keys();
    if (!keys.empty()) {
        fprintf(fd, "@<TRIPOS>PROPERTY_DATA\n");
        for (auto key : keys) {
            auto val = mol->ct(0).value(key);
            switch (val.type()) {
                case IntType:
                    fprintf(fd, "%s %lld\n", key.data(), (long long)val.asInt());
                    break;
                case FloatType:
                    fprintf(fd, "%s %f\n", key.data(), val.asFloat());
                    break;
                default:
                case StringType:
                    fprintf(fd, "%s %s\n", key.data(), val.asString().data());
                    break;
            }
        }
        fprintf(fd, "\n");
    }
}

void desres::msys::ExportMol2( SystemPtr mol, std::string const& path,
                               Provenance const& provenance,
                               IdList const& ids,
                               unsigned flags) {

    AnnotatedSystem asys(mol, AnnotatedSystem::AllowBadCharges);
    const char* mode = flags & Mol2Export::Append ? "ab" : "wb";
    FILE* fd = fopen(path.c_str(), mode);
    if (!fd) MSYS_FAIL("Could not open '" << "' for writing.");
    std::shared_ptr<FILE> dtor(fd, fclose);
    export_mol2(mol,asys,ids,fd,provenance);
}

