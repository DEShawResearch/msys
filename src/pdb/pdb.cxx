/* @COPYRIGHT@ */

#include "../pdb.hxx"
#include "../elements.hxx"
#include "../analyze.hxx"
#include "../append.hxx"
#include "../import.hxx"
#include "readpdb.h"

#include <vector>
#include <boost/algorithm/string/trim.hpp>
#include <string>
#include <math.h>
#include <errno.h>
#include "../io.hxx"

using namespace desres::msys;

static const char PDB_SPACE_GROUP[] = "pdb_space_group";
static const char PDB_Z_VALUE[] = "pdb_z_value";

static void strip_whitespace(char *buf) {
    char *ptr = buf;
    if (!buf) return;
    while (isspace(*ptr)) ++ptr;
    while (*ptr && !isspace(*ptr)) {
        *buf++ = *ptr++;
    }
    *buf='\0';
}

static void strip_nonalpha(char *buf) {
    char *ptr = buf;
    if (!buf) return;
    while (!isalpha(*ptr)) ++ptr;
    while (*ptr && isalpha(*ptr)) {
        *buf++ = *ptr++;
    }
    *buf='\0';
}

extern "C"
char* desres_msys_import_webpdb(const char* code);

static void Fclose(FILE* fp) { if (fp) fclose(fp); }

namespace {
    class iterator : public LoadIterator {
        std::shared_ptr<FILE> fd;
    public:
        explicit iterator(std::string const& path) {
            fd.reset(fopen(path.c_str(), "r"), Fclose);
            if (!fd) MSYS_FAIL("Failed opening pdb file at " << path << ": " << strerror(errno));
        }
        SystemPtr next();
    };
}

SystemPtr iterator::next() {
    char pdbstr[PDB_BUFFER_LENGTH];

    SystemPtr mol = System::create();
    mol->addCt();
    SystemImporter imp(mol);
    Id occup_id = mol->addAtomProp("occupancy", FloatType);
    Id bfactor_id = mol->addAtomProp("bfactor", FloatType);

    int indx=PDB_EOF;
    int serial, resid;
    char name[32], resname[32], chainname[32], segid[32], residstr[32];
    char insertion[4], altloc[4], element[4];
    chainname[0]=0;
    segid[0]=0;
    do {
        indx = desres_msys_read_pdb_record(fd.get(), pdbstr);

        double x, y, z, occup, beta;
        int8_t formal_charge;
        if (indx == PDB_ATOM) {

            desres_msys_get_pdb_fields(pdbstr, PDB_BUFFER_LENGTH, &serial,
                    name, resname, chainname,
                    segid, residstr, insertion, altloc, element,
                    &x, &y, &z, &occup, &beta, &formal_charge);
            strip_whitespace(name);
            strip_whitespace(resname);
            strip_whitespace(chainname);
            strip_whitespace(segid);
            strip_whitespace(element);
            strip_whitespace(altloc);
            resid = atoi(residstr);

            Id atm = imp.addAtom(chainname, segid, resid, resname, 
                    name, insertion);
            atom_t& atom = mol->atom(atm);
            atom.x = x;
            atom.y = y;
            atom.z = z;
            atom.formal_charge = formal_charge;
            if (*element) {
                atom.atomic_number = ElementForAbbreviation(element);
            } else if (*name) {
                /* guess atomic number from name.  Strip leading and
                 * trailing non-alphanumeric characters.  If the atom name
                 * and residue name match, use the atom name as the putative
                 * element symbol.  If no match is found, use the first 
                 * character.  Last resort, use the whole name. */
                strip_nonalpha(name);
                if (!strcmp(name, resname)) {
                    atom.atomic_number = ElementForAbbreviation(name);
                }
                if (atom.atomic_number==0) {
                    char tmp[2] = {name[0], 0};
                    atom.atomic_number = ElementForAbbreviation(tmp);
                }
                if (atom.atomic_number==0) {
                    atom.atomic_number = ElementForAbbreviation(name);
                }
            }
            if (strlen(altloc)) {
                Id altlocid = mol->addAtomProp("altloc", StringType);
                mol->atomPropValue(atm,altlocid) = altloc;
            }
            mol->atomPropValue(atm,occup_id) = occup;
            mol->atomPropValue(atm,bfactor_id) = beta;

        } else if (indx==PDB_CRYST1) {
            double alpha, beta, gamma, a, b, c;
            char space[12];
            int zvalue;
            desres_msys_get_pdb_cryst1(pdbstr,&alpha,&beta,&gamma,&a,&b,&c,
                                       space, &zvalue);
            space[11]='\0';
            std::string s(space);
            boost::trim(s);
            if (!s.empty()) {
                Id id = mol->ct(0).add(PDB_SPACE_GROUP, StringType);
                mol->ct(0).value(id) = s;
            }
            if (zvalue) {
                Id id= mol->ct(0).add(PDB_Z_VALUE, IntType);
                mol->ct(0).value(id) = zvalue;
            }
            ImportPDBUnitCell(a,b,c,alpha,beta,gamma,mol->global_cell[0]);

        } else if (indx==PDB_TER) {
            imp.terminateChain(chainname, segid);

        } else if (indx==PDB_END) {
            break;
        }

    } while (indx != PDB_EOF);

    if (indx==PDB_EOF && mol->maxAtomId()==0) return SystemPtr();

    GuessBondConnectivity(mol);
    Analyze(mol);

    return mol;
}

std::string desres::msys::FetchPDB(std::string const& code) {
    char* buf = desres_msys_import_webpdb(code.c_str());
    if (!buf) {
        MSYS_FAIL("Could not read pdb code " << code);
    }
    std::string s(buf);
    free(buf);
    return s;
}

SystemPtr desres::msys::ImportWebPDB(std::string const& code) {
    char* buf = desres_msys_import_webpdb(code.c_str());
    if (!buf) {
        MSYS_FAIL("Could not read pdb code " << code);
    }
    std::shared_ptr<char> ptr(buf, free);
    char temp[] = "/tmp/msys_webpdb_XXXXXX";
    int fd = mkstemp(temp);
    if (fd<0) {
        MSYS_FAIL("Could not open temporary file for webpdb");
    }
    ssize_t n = strlen(buf);
    if (::write(fd, buf, n)!=n) {
        unlink(temp);
        close(fd);
        MSYS_FAIL("Failed writing temporary file for webpdb");
    }
    ::lseek(fd,0,SEEK_SET);
    SystemPtr mol;
    try {
        mol = ImportPDB(temp);
    }
    catch (std::exception& e) {
        unlink(temp);
        close(fd);
        MSYS_FAIL("Error reading contents of webpdb " << code << " : " << e.what());
    }
    unlink(temp);
    close(fd);
    mol->name = code;
    return mol;
}

SystemPtr desres::msys::ImportPDB( std::string const& path ) {
    SystemPtr ct, mol = System::create();
    iterator it(path);
    while ((ct=it.next())) AppendSystem(mol, ct);
    mol->name = path;
    return mol;
}

LoadIteratorPtr desres::msys::PDBIterator(std::string const& path) {
    return LoadIteratorPtr(new iterator(path));
}

void desres::msys::ImportPDBCoordinates(SystemPtr mol, std::string const& path) {
    char pdbstr[PDB_BUFFER_LENGTH];
    FILE* fd = fopen(path.c_str(), "r");
    if (!fd) MSYS_FAIL("Failed opening pdb file for reading at " << path);

    int indx=PDB_EOF;
    Id id=0;
    do {
        indx = desres_msys_read_pdb_record(fd, pdbstr);
        if (indx == PDB_ATOM) {
            atom_t& atm = mol->atom(id++);
            desres_msys_get_pdb_coordinates(pdbstr, &atm.x, &atm.y, &atm.z,
                    NULL, NULL, &atm.formal_charge);
        } else if (indx==PDB_CRYST1) {
            double alpha, beta, gamma, a, b, c;
            desres_msys_get_pdb_cryst1(pdbstr,&alpha,&beta,&gamma,&a,&b,&c);
            ImportPDBUnitCell(a,b,c,alpha,beta,gamma,mol->global_cell[0]);
        }
    } while (indx!=PDB_EOF);
}

static double dotprod(const double* x, const double* y) {
    return x[0]*y[0] + x[1]*y[1] + x[2]*y[2];
}

void desres::msys::ImportPDBUnitCell(
        double a, double b, double c,
        double alpha, double beta, double gamma,
        double* m_box) {

    double A[3], B[3], C[3];

    // Convert VMD's unit cell information
    double cosBC = sin( ((90 - alpha ) / 180) * M_PI );
    double cosAC = sin( ((90 - beta  ) / 180) * M_PI );
    double cosAB = sin( ((90 - gamma ) / 180) * M_PI );
    double sinAB = cos( ((90 - gamma ) / 180) * M_PI );

    double Ax = a;
    double Ay = 0;
    double Az = 0;
    double Bx = b * cosAB;
    double By = b * sinAB;
    double Bz = 0;
    double Cx,Cy,Cz;
    if (sinAB != 0) {
        Cx = cosAC;
        Cy = (cosBC - cosAC*cosAB) / sinAB;
        Cz = sqrt(1-Cx*Cx-Cy*Cy);
        Cx *= c;
        Cy *= c;
        Cz *= c;
    } else {
        Cx=Cy=Cz=0;
    }
    A[0] = Ax; A[1] = Ay; A[2] = Az;
    B[0] = Bx; B[1] = By; B[2] = Bz;
    C[0] = Cx; C[1] = Cy; C[2] = Cz;

    /* put vectors in rows of box */
    m_box[0] = A[0]; m_box[1] = A[1]; m_box[2] = A[2];
    m_box[3] = B[0]; m_box[4] = B[1]; m_box[5] = B[2];
    m_box[6] = C[0]; m_box[7] = C[1]; m_box[8] = C[2];
}

void desres::msys::ExportPDBUnitCell(
        const double* m_box,
        double *a, double *b, double *c,
        double *alpha, double *beta, double *gamma) {

    double A[3] = { m_box[0], m_box[1], m_box[2] };
    double B[3] = { m_box[3], m_box[4], m_box[5] };
    double C[3] = { m_box[6], m_box[7], m_box[8] };

    *a = *b = *c = 1;
    *alpha = *beta = *gamma = 90;

    /* store lengths */
    *a = sqrt(dotprod(A,A));
    *b = sqrt(dotprod(B,B));
    *c = sqrt(dotprod(C,C));

    if (*a && *b && *c) {
        /* compute angles */
        double cosAB = dotprod(A,B)/(a[0] * b[0]);
        double cosAC = dotprod(A,C)/(a[0] * c[0]);
        double cosBC = dotprod(B,C)/(b[0] * c[0]);

        // clamp
        if (cosAB > 1.0) cosAB = 1.0; else if (cosAB < -1.0) cosAB = -1.0;
        if (cosAC > 1.0) cosAC = 1.0; else if (cosAC < -1.0) cosAC = -1.0;
        if (cosBC > 1.0) cosBC = 1.0; else if (cosBC < -1.0) cosBC = -1.0;

        /* convert to angles using asin to avoid nasty rounding when we are */
        /* close to 90 degree angles. */
        *alpha = 90.0 - asin(cosBC) * 90.0 / M_PI_2; /* cosBC */
        *beta  = 90.0 - asin(cosAC) * 90.0 / M_PI_2; /* cosAC */
        *gamma = 90.0 - asin(cosAB) * 90.0 / M_PI_2; /* cosAB */
    }
}

void desres::msys::ExportPDB(SystemPtr mol, std::string const& path,
                             unsigned flags) {
    FILE* fd = fopen(path.c_str(), flags & PDBExport::Append ? "ab" : "wb");
    if (!fd) MSYS_FAIL("Failed opening pdb file for writing at " << path);
    std::shared_ptr<FILE> defer_close(fd, fclose);

    if (flags & PDBExport::Append) {
        fprintf(fd, "ENDMDL\n");
    }
    
    double box[9];
    for (int i=0; i<3; i++) for (int j=0; j<3; j++) 
        box[3*i+j]=mol->global_cell[i][j];
    double A,B,C,alpha,beta,gamma;
    ExportPDBUnitCell(box, &A,&B,&C,&alpha,&beta,&gamma);
    std::string space = "P 1";
    int z = 1;
    if (mol->ctCount()>0) {
        if (mol->ct(0).has(PDB_SPACE_GROUP)) { 
            space = mol->ct(0).value(PDB_SPACE_GROUP).asString();
        }
        if (mol->ct(0).has(PDB_Z_VALUE)) { 
            z = mol->ct(0).value(PDB_Z_VALUE);
        }
    }
    fprintf(fd, "%6s%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f %-11s%4d\n",
            "CRYST1", A,B,C,alpha,beta,gamma,space.c_str(), z);

    int index=0;
    const char* chain = NULL;
    const char* resname = NULL;
    int resid = 0;
    Id bfactor_index = mol->atomPropIndex("bfactor");
    Id occupancy_index = mol->atomPropIndex("occupancy");
    for (Id chn=0; chn<mol->maxChainId(); chn++) {
        if (!mol->hasChain(chn)) continue;
        const char* segid = mol->chain(chn).segid.c_str();
        chain = mol->chain(chn).name.c_str();
        if (chn>0 && mol->chain(chn).name == mol->chain(chn-1).name) {
            desres_msys_write_ter_record(fd, ++index, resname, chain, resid);
        }
        for (Id res : mol->residuesForChain(chn)) {
            resid = mol->residue(res).resid;
            resname = mol->residue(res).name.c_str();
            const char* insertion = mol->residue(res).insertion.c_str();

            for (Id atm : mol->atomsForResidue(res)) {
                int anum = mol->atom(atm).atomic_number;
                const char* name = mol->atom(atm).name.c_str();
                const char* elementsym = AbbreviationForElement(anum);
                const char* altloc = " ";
                double x = mol->atom(atm).x;
                double y = mol->atom(atm).y;
                double z = mol->atom(atm).z;
                double occ = bad(occupancy_index) ? 1.0 :
                                mol->atomPropValue(atm, occupancy_index);
                double beta = bad(bfactor_index) ? 0.0 :
                                mol->atomPropValue(atm, bfactor_index);
                int formal_charge = mol->atom(atm).formal_charge;

                ++index;
                if (!desres_msys_write_raw_pdb_record(
                        fd, "ATOM", index, name, resname, resid,
                        insertion, altloc, elementsym, x, y, z,
                        occ, beta, formal_charge, chain, segid)) {
                    MSYS_FAIL("Failed writing PDB to " << path << " at line " << index);
                }
            }
        }
    }
}

