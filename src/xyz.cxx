#include "xyz.hxx"
#include "elements.hxx"
#include "analyze.hxx"
#include <errno.h>
#include <stdio.h>
#include <boost/algorithm/string.hpp>
#include <string.h>

namespace desres { namespace msys {

    SystemPtr ImportXYZ( std::string const& path ) {
        FILE* fd = fopen(path.c_str(), "rb");
        if (!fd) {
            MSYS_FAIL("Error opening '" << path << "' for reading: " << strerror(errno));
        }
        boost::shared_ptr<FILE> dtor(fd, fclose);
        int natoms;
        char buf[256];

        /* read #atoms */
        if (!fgets(buf, sizeof(buf), fd) || sscanf(buf, "%d", &natoms)!=1) {
            MSYS_FAIL("Failed reading number of atoms in xyz file " << path);
        }

        SystemPtr mol = System::create();
        mol->addChain();
        mol->addResidue(0);
        mol->residue(0).resid = 1;
        for (int i=0; i<natoms; i++) mol->addAtom(0);

        /* read molecule name */
        if (!fgets(buf, sizeof(buf), fd)) {
            MSYS_FAIL("Failed reading molecule name in xyz file " << path);
        }
        mol->name = buf;
        boost::trim(mol->name);

        /* read atoms */
        for (int i=0; i<natoms; i++) {
            if (!fgets(buf, sizeof(buf), fd)) {
                MSYS_FAIL("Failed reading atom " << i+1 << " in xyz file " << path);
            }
            double x,y,z;
            char name[256];
            if (4!=sscanf(buf, "%s %lf %lf %lf", name, &x, &y, &z)) {
                MSYS_FAIL("Failed parsing atom " << i+1 << " in xyz file " << path);
            }
            atom_t& atom = mol->atom(i);
            atom.x = x;
            atom.y = y;
            atom.z = z;

            /* guess atomic number: if name is a number, assume it's the
             * atomic number; otherwise treat as element name. */
            if (isalpha(name[0])) {
                atom.atomic_number = ElementForAbbreviation(name);
                atom.name = name;
            } else {
                atom.atomic_number = atoi(name);
                atom.name = AbbreviationForElement(atom.atomic_number);
            }
        }

        GuessBondConnectivity(mol);
        mol->analyze();
        return mol;
    }

}}

