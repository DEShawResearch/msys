#include "../amber.hxx"
#include <fstream>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>

using namespace desres::msys;

void desres::msys::ImportCrdCoordinates( SystemPtr mol, 
                                         std::string const& path ) {
    std::ifstream in(path);
    if (!in) MSYS_FAIL("Could not open crd file at " << path);
    std::string line;
    
    /* read the title */
    if (!std::getline(in, line)) 
        MSYS_FAIL("Reading title of crd file at " << path);
    
    /* read the number of atoms */
    if (!std::getline(in, line))
        MSYS_FAIL("Reading number of atoms of crd file at " << path);
    
    Id natoms=BadId;
    try {
        boost::trim(line);
        natoms = boost::lexical_cast<Id>(line);
    }
    catch (boost::bad_lexical_cast& e) {
        MSYS_FAIL("Could not parse number of atoms '" << line << "' in crd file " << path << ": " << e.what());
    }

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
        for (int j=0; j<2; j++) {
            atom_t& atm = mol->atom(ids[i]);
            try { 
                for (int k=0; k<3; k++) {
                    std::string s = line.substr(width*(3*j+k),width);
                    boost::trim(s);
                    (&atm.x)[k] = boost::lexical_cast<Float>(s);
                }
            }
            catch (std::exception& e) {
                MSYS_FAIL("Parsing crd file at " << path << ": " << line);
            }
            if (++i == natoms) break;
        }
    }
}    
   
