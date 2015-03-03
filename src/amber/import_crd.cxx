#include "../amber.hxx"
#include "../molfile/molfile_homebox.h"
#include <fstream>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>

using namespace desres::msys;

void desres::msys::ImportCrdCoordinates( SystemPtr mol, 
                                         std::string const& path ) {
    std::ifstream in(path.c_str());
    if (!in) MSYS_FAIL("Could not open crd file at " << path);
    std::string line;
    
    /* read the title */
    if (!std::getline(in, line)) 
        MSYS_FAIL("Reading title of crd file at " << path);
    
    /* read the number of atoms.  Use atoi so that we can read rst7 too */
    if (!std::getline(in, line))
        MSYS_FAIL("Reading number of atoms of crd file at " << path);
    Id natoms=atoi(line.c_str());

    /* number must match molecule */
    if (natoms != mol->atomCount())
        MSYS_FAIL("Number of atoms " << natoms << " in crd file " << path << " does not match number in system " << mol->atomCount());
    
    /* looks to me like we've allotted 12 columns per float */
    static const int width = 12;
    Id i=0;
    IdList ids = mol->atoms();
    while (i<natoms) {
        if (!std::getline(in, line))
            MSYS_FAIL("Reading coordinates of crd file at " << path);
        for (int j=0; j<2; j++) { /* at most two positions per line */
            atom_t& atm = mol->atom(ids[i]);
            try { 
                for (int k=0; k<3; k++) {
                    std::string s = line.substr(width*(3*j+k),width);
                    boost::trim(s);
                    (&atm.x)[k] = boost::lexical_cast<Float>(s);
                }
            }
            catch (std::exception& e) {
                MSYS_FAIL("Parsing positions in crd file at " << path << ": " << line
                        << "\n" << e.what());
            }
            if (++i == natoms) break;
        }
    }
    /* velocities */
    static const Float SCALE = 20.455;
    i=0;
    while (i<natoms) {
        if (!std::getline(in, line)) {
            if (i==0) break;
            MSYS_FAIL("Reading velocities of crd file at " << path);
        }
        for (int j=0; j<2; j++) { /* at most two positions per line */
            atom_t& atm = mol->atom(ids[i]);
            try { 
                for (int k=0; k<3; k++) {
                    std::string s = line.substr(width*(3*j+k),width);
                    boost::trim(s);
                    (&atm.vx)[k] = SCALE * boost::lexical_cast<Float>(s);
                }
            }
            catch (std::exception& e) {
                MSYS_FAIL("Parsing velocities in crd file at " << path << ": " << line
                        << "\n" << e.what());
            }
            if (++i == natoms) break;
        }
    }
    /* unit cell */
    if (std::getline(in, line)) {
        /* read a, b, c, alpha, beta, gamma */
        Float pdbbox[6];
        for (int k=0; k<6; k++) {
            std::string s = line.substr(width*k,width);
            boost::trim(s);
            try {
                pdbbox[k] = boost::lexical_cast<Float>(s);
            }
            catch (std::exception& e) {
                MSYS_FAIL("Parsing unit cell in crd file at " << path << ": " << line << "\n" << e.what());
            }
        }
        molfile_unitcell_from_pdb(mol->global_cell[0],
                pdbbox[0], pdbbox[1], pdbbox[2],
                pdbbox[3], pdbbox[4], pdbbox[5]);
    }

}    

