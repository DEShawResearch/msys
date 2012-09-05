/* @COPYRIGHT@ */

#include "../pdb.hxx"
#include "../elements.hxx"
#include "readpdb.h"

#include <vector>
#include <boost/foreach.hpp>
#include <string>

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

SystemPtr desres::msys::ImportPDB( std::string const& path ) {

    char pdbstr[PDB_BUFFER_LENGTH];

    FILE* fd = fopen(path.c_str(), "r");
    if (!fd) MSYS_FAIL("Failed opening pdb file for reading at " << path);
    boost::shared_ptr<FILE> defer_close(fd, fclose);

    SystemPtr mol = System::create();
    SystemImporter imp(mol);

    int indx=PDB_EOF;
    do {
        indx = desres_msys_read_pdb_record(fd, pdbstr);
        if (indx == PDB_ATOM) {
            int serial, resid;
            char name[32], resname[32], chainname[32], segid[32], residstr[32];
            char insertion[4], altloc[4], element[4];
            double x, y, z, occup, beta;

            desres_msys_get_pdb_fields(pdbstr, PDB_BUFFER_LENGTH, &serial,
                    name, resname, chainname,
                    segid, residstr, insertion, altloc, element,
                    &x, &y, &z, &occup, &beta);
            strip_whitespace(name);
            strip_whitespace(resname);
            strip_whitespace(chainname);
            strip_whitespace(segid);
            strip_whitespace(element);
            resid = atoi(residstr);

            Id atm = imp.addAtom(chainname, segid, resid, resname, name);
            atom_t& atom = mol->atom(atm);
            atom.x = x;
            atom.y = y;
            atom.z = z;
            atom.atomic_number = ElementForAbbreviation(element);
        }

    } while (indx != PDB_END && indx != PDB_EOF);

    mol->analyze();
    return mol;
}

void desres::msys::ExportPDB(SystemPtr mol, std::string const& path) {
    FILE* fd = fopen(path.c_str(), "w");
    if (!fd) MSYS_FAIL("Failed opening pdb file for writing at " << path);
    boost::shared_ptr<FILE> defer_close(fd, fclose);
    
    int index=0;
    for (Id chn=0; chn<mol->maxChainId(); chn++) {
        if (!mol->hasChain(chn)) continue;
        const char* segid = mol->chain(chn).segid.c_str();
        const char* chain = mol->chain(chn).name.c_str();

        BOOST_FOREACH(Id res, mol->residuesForChain(chn)) {
            int resid = mol->residue(res).resid;
            const char* resname = mol->residue(res).name.c_str();

            BOOST_FOREACH(Id atm, mol->atomsForResidue(res)) {
                const char* name = mol->atom(atm).name.c_str();
                const char* elementsym = " ";
                const char* insertion = " ";
                const char* altloc = " ";
                double x = mol->atom(atm).x;
                double y = mol->atom(atm).y;
                double z = mol->atom(atm).z;
                double occ = 1;
                double beta = 0;

                ++index;
                if (!desres_msys_write_raw_pdb_record(
                        fd, "ATOM", index, name, resname, resid,
                        insertion, altloc, elementsym, x, y, z,
                        occ, beta, chain, segid)) {
                    MSYS_FAIL("Failed writing PDB to " << path << " at line " << index);
                }
            }
        }
    }
}

