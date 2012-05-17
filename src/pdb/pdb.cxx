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

    typedef std::map<std::string,Id> ResMap;
    ResMap resmap;
    Id chn = BadId;
    Id res = BadId;

    FILE* fd = fopen(path.c_str(), "r");
    if (!fd) throw std::runtime_error("Failed opening pdb file");
    boost::shared_ptr<FILE> defer_close(fd, fclose);

    SystemPtr mol = System::create();

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

            /* Start a new chain if the chain name changes */
            if (bad(chn) || mol->chain(chn).name!=chainname) {
                chn = mol->addChain();
                mol->chain(chn).name = chainname;
                res = BadId;
            }
            /* Start a new residue if the resname or resid changes */
            if (bad(res) || mol->residue(res).resid!=resid
                    || mol->residue(res).name!=resname) {
                /* First see if this residue is in our hash */
                char key[64];
                sprintf(key, "%s - %s - %d", chainname, resname, resid);
                ResMap::const_iterator it=resmap.find(key);
                if (it==resmap.end()) {
                    res = mol->addResidue(chn);
                    mol->residue(res).resid = resid;
                    mol->residue(res).name = resname;
                    resmap[key]=res;
                }
            }
            /* Append particles to the end of the residue */
            Id atm = mol->addAtom(res);
            mol->atom(atm).name = name;
            mol->atom(atm).x = x;
            mol->atom(atm).y = y;
            mol->atom(atm).z = z;
        }

    } while (indx != PDB_END && indx != PDB_EOF);

    return mol;
}

