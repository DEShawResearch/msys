/* @COPYRIGHT@ */

#include "../pdb.hxx"
#include "../elements.hxx"
#include "../analyze.hxx"
#include "readpdb.h"

#include <vector>
#include <boost/foreach.hpp>
#include <string>
#include <math.h>

using namespace desres::msys;

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

SystemPtr desres::msys::ImportPDB( std::string const& path ) {

    char pdbstr[PDB_BUFFER_LENGTH];

    FILE* fd = fopen(path.c_str(), "r");
    if (!fd) MSYS_FAIL("Failed opening pdb file for reading at " << path);
    boost::shared_ptr<FILE> defer_close(fd, fclose);

    SystemPtr mol = System::create();
    SystemImporter imp(mol);
    Id occup_id = mol->addAtomProp("occupancy", FloatType);
    Id bfactor_id = mol->addAtomProp("bfactor", FloatType);

    int indx=PDB_EOF;
    int ct=0;
    int ctcount=0;
    do {
        indx = desres_msys_read_pdb_record(fd, pdbstr);
        if (indx == PDB_ATOM) {
            int serial, resid;
            char name[32], resname[32], chainname[32], segid[32], residstr[32];
            char insertion[4], altloc[4], element[4];
            double x, y, z, occup, beta;
            int formal_charge;

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
                    name, insertion, ct);
            ++ctcount;
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
            desres_msys_get_pdb_cryst1(pdbstr,&alpha,&beta,&gamma,&a,&b,&c);
            ImportPDBUnitCell(a,b,c,alpha,beta,gamma,mol->global_cell[0]);

        } else if (/* indx==PDB_TER || */ indx==PDB_END) {
            if (ctcount>0) {
                ++ct;
                ctcount=0;
            }
        }

    } while (indx != PDB_EOF);

    GuessBondConnectivity(mol);
    mol->analyze();

    return mol;
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

void desres::msys::ExportPDB(SystemPtr mol, std::string const& path) {
    FILE* fd = fopen(path.c_str(), "w");
    if (!fd) MSYS_FAIL("Failed opening pdb file for writing at " << path);
    boost::shared_ptr<FILE> defer_close(fd, fclose);
    
    double box[9];
    for (int i=0; i<3; i++) for (int j=0; j<3; j++) 
        box[3*i+j]=mol->global_cell[i][j];
    double A,B,C,alpha,beta,gamma;
    ExportPDBUnitCell(box, &A,&B,&C,&alpha,&beta,&gamma);

    fprintf(fd, "%6s%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f %-11s%4d\n",
            "CRYST1", A,B,C,alpha,beta,gamma,"P 1", 1);

    int index=0;
    for (Id chn=0; chn<mol->maxChainId(); chn++) {
        if (!mol->hasChain(chn)) continue;
        const char* segid = mol->chain(chn).segid.c_str();
        const char* chain = mol->chain(chn).name.c_str();

        BOOST_FOREACH(Id res, mol->residuesForChain(chn)) {
            int resid = mol->residue(res).resid;
            const char* resname = mol->residue(res).name.c_str();

            BOOST_FOREACH(Id atm, mol->atomsForResidue(res)) {
                int anum = mol->atom(atm).atomic_number;
                const char* name = mol->atom(atm).name.c_str();
                const char* elementsym = AbbreviationForElement(anum);
                const char* insertion = " ";
                const char* altloc = " ";
                double x = mol->atom(atm).x;
                double y = mol->atom(atm).y;
                double z = mol->atom(atm).z;
                double occ = 1;
                double beta = 0;
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

