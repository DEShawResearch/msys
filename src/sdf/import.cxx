#include "../sdf.hxx"
#include "../sssr.hxx"
#include "elements.hxx"
#include <boost/foreach.hpp>
#include <boost/algorithm/string.hpp> /* for boost::trim */
#include <boost/circular_buffer.hpp>
#include <stdio.h>
#include <errno.h>
#include <string>

#include <fstream>


using namespace desres::msys;

namespace {
    int stringToInt(std::string const& str){
        char* stop;
        int res = strtol( str.c_str(), &stop, 10 );
        if ( *stop != 0 ) MSYS_FAIL("Bad int Specification:\n" << str);
        // printf("stringToInt: '%s' -> %d\n",str.c_str(),res);
        return res;
    }
    
    double stringToDouble(std::string const& str){
        char* stop;
        double res = strtod( str.c_str(), &stop );
        if ( *stop != 0 ) MSYS_FAIL("Bad double Specification:\n" << str);
        // printf("stringToDouble: '%s' -> %f\n",str.c_str(),res);
        return res;
    }

    SystemPtr parseV2000mol(std::istream &ifs, 
                            std::string const& name, 
                            std::string const& program, 
                            std::string const& comment,
                            int natoms, int nbonds){
        /* Create Molecule */
        SystemPtr mol = System::create();
        Id chn = mol->addChain();
        Id res = mol->addResidue(chn);
        // mol->residue(res).name = resname;
        mol->name = name;

        IdList atoms;
        std::vector<int> valenceList;

        std::string line;
        /* Load Atoms */
        for (int i=0; i<natoms; i++) {
            if (!std::getline(ifs, line)) {
                MSYS_FAIL("Missing expected Atom record " << i+1);
            }

            /* Allow atom lines to have incomplete info (must have up to element) */
            if(line.size()<34){
                MSYS_FAIL("Malformed atom line:\n"<<line);
            }
        
            Id atm = mol->addAtom(res);
            atoms.push_back(atm);
            atom_t& atom = mol->atom(atm);
        
            atom.x = stringToDouble(line.substr(0,10));
            atom.y = stringToDouble(line.substr(10,10));
            atom.z = stringToDouble(line.substr(20,10));
            std::string sym=line.substr(31,3);
            boost::trim(sym);
            atom.name = sym;
            atom.atomic_number = ElementForAbbreviation(sym.c_str());

            if(line.size()>=39){
                int q=stringToInt(line.substr(36,3));
                /* Values of parsed charge field: 
                   0 = uncharged or value other than these:
                   1 = +3, 2 = +2, 3 = +1,
                   4 = doublet radical, 
                   5 = -1, 6 = -2, 7 = -3 
                */
                if(q<0 || q>7){
                    MSYS_FAIL("charge specification for atom is out of range:\n" << line);
                }
                if(q==4){
                    fprintf(stderr, "WARNING: Treating doublet radical charge specification as q=0\n");
                }
                atom.charge = q==0 ? 0 : 4-q;
                
                if(line.size()>=51){
                    int v=stringToInt(line.substr(48,3));
                    /* Values of parsed valence field: 
                       0 = no marking (default)
                       (1 to 14) = (1 to 14) 
                       15 = zero valence
                    */
                    if(v<0 || v>15){
                        MSYS_FAIL("valence specification for atom is out of range:\n" << line);
                    }
                    /* Use valence list to check bonds? */
                    valenceList.push_back(v);
                }
            }
        } 

        /* Load Bonds */
        for (int i=0; i<nbonds; i++) {

            if (!std::getline(ifs, line)) {
                MSYS_FAIL("Missing expected Bond record " << i+1);
            }
            if(line.size()<9){
                MSYS_FAIL("Malformed bond line:\n"<<line);
            }

            int ai=stringToInt(line.substr(0,3));
            int aj=stringToInt(line.substr(3,3));
            int type=stringToInt(line.substr(6,3));

            if (ai<1 || aj<1 || ai>natoms || aj>natoms) {
                MSYS_FAIL("Invalid atom reference in Bond record:\n" << line);
            }
            if(type<1 || type>4){
                MSYS_FAIL("Invalid bond type in Bond record:\n" << line);
            }
            Id bnd = mol->addBond(atoms.at(ai-1), atoms.at(aj-1));
            bond_t& bond = mol->bond(bnd);
            if (type==4) {
                bond.order = 1;
                bond.resonant_order = 1.5;       
            }else{
                bond.order = bond.resonant_order = type;
            }
        }
        BOOST_FOREACH(IdList const& ring, 
                      GetSSSR(mol, mol->atoms())) {
            bool all_resonant = true;
            for (Id i=1; i<ring.size(); i++) {
                Id bid = mol->findBond(ring[i-1], ring[i]);
                if (bad(bid)) MSYS_FAIL("Missing bond!");
                if (mol->bond(bid).resonant_order!=1.5) {
                    all_resonant = false;
                    break;
                }
            }
            if (all_resonant) {
                for (Id i=1; i<ring.size(); i++) {
                    Id bid = mol->findBond(ring[i-1], ring[i]);
                    mol->bond(bid).order = i % 2;
                }
            }
        }

        /* properties, data and end of molecule delimeter scanning 
           Currently not parsed */
        while (std::getline(ifs, line)) {
            if(line.compare(0,3,"M  ")==0){
                if(line.compare(3,3,"CHG")==0){
                    // Atom Charge Info
                }else if(line.compare(3,3,"END")==0){
                    // End of molecule properties block
                }
            }else if(line.compare(0,3,"A  ")==0 || line.compare(0,3,"G  ")==0){
                std::getline(ifs, line);
            }else if(line.compare(0,3,"V  ")==0){
            }else if(line.compare(0,2,"> ")==0){
                while (std::getline(ifs, line)){
                    boost::trim(line);
                    if(line=="")break;
                }
            }else if (line.compare(0,4,"$$$$")==0){
                // End of molecule
                break;
            }else{
                /* According to the spec, we shouldnt find any other blank lines,
                   but we will be forgiving just in case */
                boost::trim(line);
                if(line != "")MSYS_FAIL("Unknown line in SDF file:\n"<<line);
            }
        }
        return mol;
    }

    SystemPtr parseV3000mol(std::istream &ifs, 
                            std::string const& name, 
                            std::string const& program, 
                            std::string const& comment,
                            int natoms, int nbonds){
        fprintf(stderr, "WARNING: Skipping V3000 formatted molecule\n");
        std::string line;
        while (std::getline(ifs, line)) {
            if(line.compare(0,4,"$$$$")==0) break;
        }
        return System::create();
    }

}

SystemPtr desres::msys::ImportSdf(std::string const& path) {
    std::vector<SystemPtr> mols = ImportSdfMany(path);
    if (mols.empty()) return SystemPtr();
    return mols[0];
}

std::vector<SystemPtr>
desres::msys::ImportSdfMany(std::string const& path) {

    std::ifstream ifs(path);
    if (!ifs) MSYS_FAIL("Could not open MDL MOL/SDF file for reading at " << path);

    std::vector<SystemPtr> mols;
    /* Circular buffer for header block. 
       More forgiving of too few/many blank lines in poorly formatted files */
    boost::circular_buffer<std::string> cb(3,"");

    /* scan for the counts lines */
    std::string line;
    while(std::getline(ifs, line)){
        boost::trim_right(line);
        if(line.size()==38 || line.size()==39){
            std::string version=line.substr(33,6);
            boost::trim(version);
            
            if(version=="V2000" || version=="V3000"){
                int natoms=stringToInt(line.substr(0,3));
                int nbonds=stringToInt(line.substr(3,3));
                
                /* header fields are:
                   cb[0] == name
                   cb[1] == program information
                   cb[2] == comment
                */
                if(version=="V2000"){
                    mols.push_back(parseV2000mol(ifs, cb[0], cb[1], cb[2], natoms, nbonds));
                }else if (version=="V3000") {
                    mols.push_back(parseV3000mol(ifs, cb[0], cb[1], cb[2], natoms, nbonds));
                }
                /* reset the circular buffer after a sucessful molecule parse */
                cb.assign(3,"");
            }else{
                cb.push_back(line);
            }
        }else{
            cb.push_back(line);     
        }
    }
   
    if (!ifs.eof()) MSYS_FAIL("Error reading from " << path << ": " << strerror(errno));

    printf("Loaded %lu molecules from file %s\n",mols.size(),path.c_str());
    return mols;
}
