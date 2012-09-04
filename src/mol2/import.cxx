#include "../mol2.hxx"
#include "elements.hxx"
#include <boost/foreach.hpp>
#include <boost/algorithm/string.hpp> /* for boost::trim */
#include <stdio.h>
#include <string.h>
#include <errno.h>

using namespace desres::msys;

namespace {
    enum State { Start, Skip, Molecule, Atom, Bond, Substructure };
}

SystemPtr desres::msys::ImportMol2(std::string const& path) {
    std::vector<SystemPtr> mols = ImportMol2Many(path);
    if (mols.empty()) return SystemPtr();
    return mols[0];
}

std::vector<SystemPtr>
desres::msys::ImportMol2Many(std::string const& path) {
    FILE* fd = fopen(path.c_str(), "rb");
    if (!fd) MSYS_FAIL("Could not open mol2 file for reading at " << path);
    boost::shared_ptr<FILE> dtor(fd, fclose);

    std::vector<SystemPtr> mols;
    char buf[256];
    State state = Start;
    IdList atoms;
    int natoms, nbonds, nsub;
    SystemPtr mol;
    Id chn = BadId;
    while (fgets(buf, sizeof(buf), fd)) {
        if (buf[0]=='@') {
            if (!strncmp(buf, "@<TRIPOS>MOLECULE", 17)) {
                state = Molecule;
            } else if (!strncmp(buf, "@<TRIPOS>ATOM", 13)) {
                state = Atom;
            } else if (!strncmp(buf, "@<TRIPOS>BOND", 13)) {
                state = Bond;
            } else if (!strncmp(buf, "@<TRIPOS>SUBSTRUCTURE", 21)) {
                state = Substructure;
            } else {
                state = Start;
            }
        }
        switch (state) {
            case Molecule: 
                mol = System::create();
                mols.push_back(mol);
                chn = mol->addChain();
                atoms.clear();
                /* read mol_name */
                fgets(buf, sizeof(buf), fd); 
                mol->name = buf;
                boost::trim(mol->name);
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
                /* skip the rest */
                state = Skip;
                break;

            case Skip: break;

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
                    if (rc>=7 && (resid<1)) {
                        fprintf(stderr, "WARNING: Invalid subst_id %d will be assumed to be 1 in:\n%s\n",
                                resid, buf);
                        resid=1;
                    }
                    int nres = mol->residueCountForChain(chn);
                    for (; nres<resid; ++nres) {
                        Id res = mol->addResidue(chn);
                        mol->residue(res).resid = nres+1;
                        mol->residue(nres).name = resname;
                    }
                    Id res = mol->residuesForChain(chn).at(resid-1);
                    Id atm = mol->addAtom(res);
                    atoms.push_back(atm);
                    atom_t& atom = mol->atom(atm);
                    atom.name = name;
                    atom.x = x;
                    atom.y = y;
                    atom.z = z;
                    atom.charge = q;
                    char* dot = strchr(type, '.');
                    if (dot) *dot='\0';
                    atom.atomic_number = ElementForAbbreviation(type);
                    mol->residue(res).name = resname;
                }
                state = Skip;
                break;

            case Bond:
                for (int i=0; i<nbonds; i++) {
                    int ai, aj;
                    char btype[32];
                    if (!fgets(buf, sizeof(buf), fd)) {
                        MSYS_FAIL("Missing expected Bond record " << i+1);
                    }
                    if (sscanf(buf, "%*d %d %d %s", &ai, &aj, btype)!=3) {
                        MSYS_FAIL("Could not parse Bond record:\n" << buf);
                    }
                    if (ai<1 || aj<1 || ai>natoms || aj>natoms) {
                        MSYS_FAIL("Invalid atoms in Bond record:\n" << buf);
                    }
                    mol->addBond(atoms.at(ai-1), atoms.at(aj-1));
                }
                state = Skip;
                break;

            case Substructure:
                for (int i=0; i<nsub; i++) {
                    if (!fgets(buf, sizeof(buf), fd)) {
                        MSYS_FAIL("Missing expected Substructure record " << i+1);
                    }

                }
            default: ;
        }
    }
    if (!feof(fd)) MSYS_FAIL("Error reading from " << path << ": " << strerror(errno));
    return mols;
}
