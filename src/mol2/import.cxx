#include "../mol2.hxx"
#include "../sssr.hxx"
#include "../import.hxx"
#include "../analyze.hxx"
#include "../elements.hxx"
#include "../append.hxx"
#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <unordered_map>

using namespace desres::msys;

namespace {
    class iterator : public LoadIterator {
        FILE* fd;
        enum State { Skip, Molecule, Atom, Bond, Substructure } state;
        char buf[256];
        long file_offset;

        /* go to next molecule */
        void advance();
    public:
        ~iterator() {
            if (fd) fclose(fd);
        }
        explicit iterator(std::string const& path) : fd(), file_offset() {
            fd = fopen(path.c_str(), "rb");
            if (!fd) {
                MSYS_FAIL("Could not open mol2 file for reading at " << path);
            }
            advance();
        }
        SystemPtr next();
    };
}

SystemPtr desres::msys::ImportMol2(std::string const& path) {
    SystemPtr ct, mol = System::create();
    iterator it(path);
    while ((ct=it.next())) AppendSystem(mol,ct);
    mol->name = path;
    return mol;
}

std::vector<SystemPtr>
desres::msys::ImportMol2Many(std::string const& path) {
    iterator iter(path);
    std::vector<SystemPtr> mols;
    SystemPtr mol;
    while ((mol=iter.next())) {
        mols.push_back(mol);
    }
    return mols;
}

LoadIteratorPtr desres::msys::Mol2Iterator(std::string const& path) {
    return LoadIteratorPtr(new iterator(path));
}

void iterator::advance() {
    while (fgets(buf, sizeof(buf), fd)) {
        if (!strncmp(buf, "@<TRIPOS>MOLECULE", 17)) {
            file_offset = ftell(fd) - strlen(buf);
            state = Molecule;
            break;
        }
    }
}

SystemPtr iterator::next() {

    /* if we hit eof before reaching a Molecule section last time,
     * then return nothing but empty systems */
    if (state!=Molecule) return SystemPtr();

    int natoms, nbonds, nsub;

    SystemPtr mol = System::create();
    mol->addCt();
    mol->ct(0).add("msys_file_offset", IntType);
    mol->ct(0).value("msys_file_offset") = file_offset;
    mol->addChain(0);
    std::unordered_map<std::string, Id> chain_to_id;

    /* read mol_name */
    if (!fgets(buf, sizeof(buf), fd)) MSYS_FAIL(strerror(errno));
    mol->name = buf;
    trim(mol->name);
    mol->ct(0).setName(mol->name);
    /* read natoms, nbonds, nsub */
    if (!fgets(buf, sizeof(buf), fd)) MSYS_FAIL(strerror(errno));
    if (sscanf(buf, "%d %d %d", &natoms, &nbonds, &nsub)!=3) {
        MSYS_FAIL("Could not parse counts from line:\n" << buf);
    }
    /* read mol_type and ignore */
    if (!fgets(buf, sizeof(buf), fd)) MSYS_FAIL(strerror(errno));
    /* read charge_type and ignore */
    if (!fgets(buf, sizeof(buf), fd)) MSYS_FAIL(strerror(errno));

    while (fgets(buf, sizeof(buf), fd)) {
        if (buf[0]=='@') {
            if (!strncmp(buf, "@<TRIPOS>MOLECULE", 17)) {
                file_offset = ftell(fd) - strlen(buf);
                state = Molecule;
                break;
            } else if (!strncmp(buf, "@<TRIPOS>ATOM", 13)) {
                state = Atom;
            } else if (!strncmp(buf, "@<TRIPOS>BOND", 13)) {
                state = Bond;
            } else if (!strncmp(buf, "@<TRIPOS>SUBSTRUCTURE", 21)) {
                state = Substructure;
            } else {
                state = Skip;
            }
        }
        switch (state) {
            case Skip: 
                break;
            case Atom:
                {
                for (int i=0; i<natoms; i++) {
                    int id, resid=0;
                    Id residue = BadId, residue_id = BadId;
                    Float x,y,z,q=0;
                    char name[32], type[32], resname[32];
                    if (!fgets(buf, sizeof(buf), fd)) {
                        MSYS_FAIL("Missing expected Atom record " << i+1);
                    }
                    int rc = sscanf(buf, "%d %s %lf %lf %lf %s %u %s %lf",
                            &id, name, &x, &y, &z, type, &residue, resname, &q);
                    if (rc<8) {
                        MSYS_FAIL("Could not parse Atom record:\n" << buf);
                    }
                    // maybe the resid is encoded as digits at the end
                    // of the resname.
                    resid = residue;
                    {
                        char* p=resname;
                        while (*p && isalpha(*p)) ++p;
                        if (*p) {
                            char* end=p;
                            int newresid = strtol(p, &end, 10);
                            if (*end=='\0') {
                                resid = newresid;
                                *p = '\0';
                            }
                        }
                    }
                    if (residue == mol->residueCount() + 1) {
                        residue_id = mol->addResidue(0);
                        auto& res = mol->residueFAST(residue_id);
                        res.name = resname;
                        res.resid = resid;
                    } else if (residue == mol->residueCount()) {
                        residue_id = residue - 1;
                    } else {
                        MSYS_FAIL("Nonconsecutive residue id: " << buf);
                    }
                    Id atmid = mol->addAtom(residue_id);
                    auto& atm = mol->atomFAST(atmid);
                    atm.name = name;
                    atm.x = x;
                    atm.y = y;
                    atm.z = z;
                    atm.charge = q;
                    char* dot = strchr(type, '.');
                    if (dot) *dot='\0';
                    atm.atomic_number = ElementForAbbreviation(type);
                }
                state = Skip;
                }
                break;

            case Bond:
                for (int i=0; i<nbonds; i++) {
                    int ai, aj;
                    char type[32];
                    if (!fgets(buf, sizeof(buf), fd)) {
                        MSYS_FAIL("Missing expected Bond record " << i+1);
                    }
                    if (sscanf(buf, "%*d %d %d %s", &ai, &aj, type)!=3) {
                        MSYS_FAIL("Could not parse Bond record:\n" << buf);
                    }
                    if (ai<1 || aj<1 || ai>natoms || aj>natoms) {
                        MSYS_FAIL("Invalid atoms in Bond record:\n" << buf);
                    }
                    Id bnd = mol->addBond(ai-1, aj-1);
                    bond_t& bond = mol->bond(bnd);
                    if (!strcmp(type, "1")) {
                        bond.order = 1;
                    } else if (!strcmp(type, "2")) {
                        bond.order = 2;
                    } else if (!strcmp(type, "3")) {
                        bond.order = 3;
                    } else if (!strcmp(type, "am")) {
                        bond.order = 1;
                    } else if (!strcmp(type, "ar")) {
                        bond.order = 1;
                        bond.aromatic = true;
                    }
                }
                state = Skip;
                break;

            case Substructure:
                for (int i=0; i<nsub; i++) {
                    if (!fgets(buf, sizeof(buf), fd)) {
                        MSYS_FAIL("Missing expected Substructure record " << i+1);
                    }
                    // check that we read a consecutive id
                    Id id = BadId;
                    char chn[32];
                    memset(chn, 0, sizeof(chn));
                    if (sscanf(buf, "%u %*s %*d %*s %*d %s", &id, chn) < 1) {
                        MSYS_FAIL("Could not parse Substructure record:\n" << buf);
                    }
                    if (id != (unsigned)(i+1)) {
                        MSYS_FAIL("Nonconsecutive substructure ids:\n" << buf);
                    }
                    if (i==0) {
                        // first residue goes in first chain
                        chain_to_id[chn] = 0;
                    } else {
                        Id chnid = BadId;
                        if (chain_to_id.count(chn) == 1) {
                            chnid = chain_to_id[chn];
                        } else {
                            // first time seeing this chain.  Make a new chain.
                            chnid = mol->addChain(0);
                            chain_to_id[chn] = chnid;
                        }
                        if (chnid != 0) {
                            // this residue got added to chain 0, but it's actually in
                            // a different chain.  Fix it up.
                            mol->setChain(id-1, chnid);
                        }
                    }
                }
                // update chain info
                for (auto& p : chain_to_id) {
                    auto& name = p.first;
                    Id chnid = p.second;
                    auto& chn = mol->chain(chnid);
                    chn.name = name;
                    for (const char* c = name.data(); *c; c++) {
                        if (isdigit(*c)) {
                            chn.segid = c;
                            chn.name = std::string(name.data(), c);
                            break;
                        }
                    }
                }
                state = Skip;
            default: ;
        }
    }
    Analyze(mol);
    return mol;
}
