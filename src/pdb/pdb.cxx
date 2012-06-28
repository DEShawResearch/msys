/* @COPYRIGHT@ */

#include "../pdb.hxx"
#include "readpdb.h"

#include <vector>
#include <string>
#include <map>
#include <stdexcept>

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
    if (!fd) throw std::runtime_error("Failed opening pdb file");
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
            float x, y, z, occup, beta;

            desres_msys_get_pdb_fields(pdbstr, PDB_BUFFER_LENGTH, &serial,
                    name, resname, chainname,
                    segid, residstr, insertion, altloc, element,
                    &x, &y, &z, &occup, &beta);
            strip_whitespace(name);
            strip_whitespace(resname);
            strip_whitespace(chainname);
            strip_whitespace(segid);
            resid = atoi(residstr);

            Id atm = imp.addAtom(chainname, segid, resid, resname, name);
            mol->atom(atm).x = x;
            mol->atom(atm).y = y;
            mol->atom(atm).z = z;
        }

    } while (indx != PDB_END && indx != PDB_EOF);

    return mol;
}

