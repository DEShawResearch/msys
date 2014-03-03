#include "../mol2.hxx"
#include "../sssr.hxx"
#include "elements.hxx"
#include "append.hxx"
#include <boost/foreach.hpp>
#include <boost/algorithm/string.hpp> /* for boost::trim */
#include <stdio.h>
#include <string.h>
#include <errno.h>

using namespace desres::msys;

namespace {
    class iterator : public LoadIterator {
        FILE* fd;
        enum State { Skip, Molecule, Atom, Bond, Substructure } state;
        char buf[256];

        /* go to next molecule */
        void advance();
    public:
        ~iterator() {
            if (fd) fclose(fd);
        }
        explicit iterator(std::string const& path) : fd() {
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
    Id atype=mol->addAtomProp("sybyl_type", StringType);
    Id btype=mol->addBondProp("sybyl_type", StringType);
    Id chn = mol->addChain();

    /* read mol_name */
    fgets(buf, sizeof(buf), fd); 
    mol->name = buf;
    boost::trim(mol->name);
    mol->ct(0).setName(mol->name);
    /* read natoms, nbonds, nsub */
    natoms = nbonds = 0;
    nsub = 0;
    fgets(buf, sizeof(buf), fd);
    if (sscanf(buf, "%d %d %d", &natoms, &nbonds, &nsub)<1) {
        MSYS_FAIL("Could not parse num_atoms from line:\n" << buf);
    }
    /* read mol_type and ignore */
    fgets(buf, sizeof(buf), fd);
    /* read charge_type and ignore */
    fgets(buf, sizeof(buf), fd);

    while (fgets(buf, sizeof(buf), fd)) {
        if (buf[0]=='@') {
            if (!strncmp(buf, "@<TRIPOS>MOLECULE", 17)) {
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
                for (int i=0; i<natoms; i++) {
                    int id, resid=0;
                    Float x,y,z,q=0;
                    char name[32], type[32], resname[32];
                    if (!fgets(buf, sizeof(buf), fd)) {
                        MSYS_FAIL("Missing expected Atom record " << i+1);
                    }
                    int rc = sscanf(buf, "%d %s %lf %lf %lf %s %d %s %lf",
                            &id, name, &x, &y, &z, type, &resid, resname, &q);
                    if (rc<6) {
                        MSYS_FAIL("Could not parse Atom record:\n" << buf);
                    }
                    int nres = mol->residueCountForChain(chn);
                    for (; nres<resid; ++nres) {
                        Id res = mol->addResidue(chn);
                        mol->residue(res).resid = nres+1;
                        mol->residue(nres).name = resname;
                        // maybe the resid is encoded as digits at the end
                        // of the resname.
                        {
                            char* p=resname;
                            while (*p && isalpha(*p)) ++p;
                            if (*p) {
                                char* end=p;
                                int resid = strtol(p, &end, 10);
                                if (*end=='\0') {
                                    mol->residue(res).resid = resid;
                                    *p = '\0';
                                    mol->residue(res).name = resname;
                                }
                            }
                        }
                    }
                    Id res = mol->residuesForChain(chn).at(resid-1);
                    Id atm = mol->addAtom(res);
                    mol->atomPropValue(atm,atype)=type;
                    atom_t& atom = mol->atom(atm);
                    atom.name = name;
                    atom.x = x;
                    atom.y = y;
                    atom.z = z;
                    atom.charge = q;
                    char* dot = strchr(type, '.');
                    if (dot) *dot='\0';
                    atom.atomic_number = ElementForAbbreviation(type);
                }
                state = Skip;
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
                    mol->bondPropValue(bnd,btype)=type;
                    bond_t& bond = mol->bond(bnd);
                    if (!strcmp(type, "1")) {
                        bond.order = bond.resonant_order = 1;
                    } else if (!strcmp(type, "2")) {
                        bond.order = bond.resonant_order = 2;
                    } else if (!strcmp(type, "3")) {
                        bond.order = bond.resonant_order = 3;
                    } else if (!strcmp(type, "am")) {
                        bond.order = bond.resonant_order = 1;
                    } else if (!strcmp(type, "ar")) {
                        bond.order = 1;
                        bond.resonant_order = 1.5;
                    }
                }
                state = Skip;
                break;

            case Substructure:
                for (int i=0; i<nsub; i++) {
                    if (!fgets(buf, sizeof(buf), fd)) {
                        MSYS_FAIL("Missing expected Substructure record " << i+1);
                    }

                }
                state = Skip;
            default: ;
        }
    }
    mol->analyze();
    return mol;
}
