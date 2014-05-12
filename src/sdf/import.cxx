#include "../sdf.hxx"
#include "elements.hxx"
#include "../append.hxx"
#include <boost/foreach.hpp>
#include <boost/algorithm/string.hpp> /* for boost::trim */
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>

#include <stdio.h>
#include <errno.h>
#include <string>
#include <math.h>

#include <fstream>
#include <sstream>


using namespace desres::msys;

namespace {
    /* FIXME: these are also implemented in atomsel.. */
    int stringToInt(std::string const& str){
        char* stop;
        int res = strtol( str.c_str(), &stop, 10 );
        if ( *stop != 0 ) MSYS_FAIL("Bad int Specification: '" << str << "'");
        return res;
    }
    
    double stringToDouble(std::string const& str){
        char* stop;
        double res = strtod( str.c_str(), &stop );
        if ( *stop != 0 ) MSYS_FAIL("Bad double Specification:\n" << str);
        return res;
    }

    void add_typed_keyval(String const& key, String const& val, 
                          component_t& ct) {
        try {
            int v = stringToInt(val);
            ct.add(key,IntType);
            ct.value(key)=v;
            return;
        } catch (Failure& e) {
        }
        try {
            double v = stringToDouble(val);
            if (isfinite(v)) {
                ct.add(key,FloatType);
                ct.value(key)=v;
                return;
            }
        } catch (Failure& e) {
        }
        ct.add(key,StringType);
        ct.value(key)=val;
    }

    class iterator : public LoadIterator {
        std::ifstream file;
        boost::iostreams::filtering_istream in;

        void skip_to_end();

    public:
        explicit iterator(std::string const& path)
        : file(path.c_str()) {
            if (!file) {
                MSYS_FAIL("Failed opening SDF file for reading at " << path);
            }
            init(file);
        }

        void init(std::istream& input) {
            /* FIXME - copied from mae/mae.cxx */
            /* check for gzip magic number */
            if (input.get()==0x1f && input.get()==0x8b) {
                in.push(boost::iostreams::gzip_decompressor());
            }
            input.seekg(0);
            in.push(input);
        }

        SystemPtr next();
    };
}

void iterator::skip_to_end() {
    std::string line;
    while (std::getline(in,line)) {
        if (line.compare(0,4,"$$$$")==0) break;
    }
}

SystemPtr iterator::next() {
    if (in.eof()) return SystemPtr();

    /* Header block contains exactly three lines.  */
    std::string h1, h2, h3;
    std::getline(in,h1);
    std::getline(in,h2);
    std::getline(in,h3);
    if (!in) {
        if (in.eof()) return SystemPtr();
        skip_to_end();
        MSYS_FAIL("Failed reading header block.");
    }

    /* create system; assign molecule name from line 1*/
    SystemPtr mol = System::create();
    Id chn = mol->addChain();
    Id res = mol->addResidue(chn);
    boost::trim(h1);
    mol->name = h1;
    mol->ct(0).setName(h1);

    /* parse counts line */
    std::string line;
    if (!std::getline(in,line)) {
        skip_to_end();
        MSYS_FAIL("Failed reading counts");
    }

    int natoms = stringToInt(line.substr(0,3));
    int nbonds = stringToInt(line.substr(3,3));

    std::vector<int> valenceList;

    /* Load Atoms */
    for (int i=0; i<natoms; i++) {
        if (!std::getline(in, line)) {
            skip_to_end();
            MSYS_FAIL("Missing expected Atom record " << i+1);
        }

        /* Allow atom lines to have incomplete info (must have up to element) */
        if(line.size()<34){
            skip_to_end();
            MSYS_FAIL("Malformed atom line:\n"<<line);
        }
    
        Id atm = mol->addAtom(res);
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
                skip_to_end();
                MSYS_FAIL("charge specification for atom is out of range:\n" << line);
            }
            if(q==4){
                skip_to_end();
                fprintf(stderr, "WARNING: Treating doublet radical charge specification as q=0\n");
            }
            /* Just set formal charge, not 'charge', for consistency with
             * the exporter, which only looks at formal charge.  */
            atom.formal_charge = q==0 ? 0 : 4-q;
            
            if(line.size()>=51){
                int v=stringToInt(line.substr(48,3));
                /* Values of parsed valence field: 
                   0 = no marking (default)
                   (1 to 14) = (1 to 14) 
                   15 = zero valence
                */
                if(v<0 || v>15){
                    skip_to_end();
                    MSYS_FAIL("valence specification for atom is out of range:\n" << line);
                }
                /* Use valence list to check bonds? */
                valenceList.push_back(v);
            }
        }
    } 

    /* Load Bonds */
    for (int i=0; i<nbonds; i++) {

        if (!std::getline(in, line)) {
            skip_to_end();
            MSYS_FAIL("Missing expected Bond record " << i+1);
        }
        if(line.size()<9){
            skip_to_end();
            MSYS_FAIL("Malformed bond line:\n"<<line);
        }

        int ai=stringToInt(line.substr(0,3));
        int aj=stringToInt(line.substr(3,3));
        int type=stringToInt(line.substr(6,3));

        if (ai<1 || aj<1 || ai>natoms || aj>natoms) {
            skip_to_end();
            MSYS_FAIL("Invalid atom reference in Bond record:\n" << line);
        }
        if(type<1 || type>4){
            skip_to_end();
            MSYS_FAIL("Invalid bond type in Bond record:\n" << line);
        }
        Id bnd = mol->addBond(ai-1, aj-1);
        bond_t& bond = mol->bond(bnd);
        if (type==4) {
            bond.order = 1;
            bond.resonant_order = 1.5;       
        }else{
            bond.order = bond.resonant_order = type;
        }
    }
    /* properties, data and end of molecule delimeter scanning.
     * Incompletely parsed. */
    bool cleared_m_chg = false;
    while (std::getline(in, line)) {
        if (line.compare(0,3,"M  ")==0){
            if (line.compare(3,3,"CHG")==0){
                /* Atom Charge Info supersedes atom record lines */
                if (!cleared_m_chg) {
                    cleared_m_chg = true;
                    for (Id i=0; i<mol->atomCount(); i++) {
                        mol->atom(i).formal_charge=0;
                    }
                }
                std::istringstream ss(line.substr(6));
                int i, n;
                ss >> n;
                for (i=0; i<n; i++) {
                    int ai,fc;
                    ss >> ai >> fc;
                    if (!ss) {
                        skip_to_end();
                        MSYS_FAIL("Failed parsing CHG line " << line);
                    }
                    mol->atom(ai-1).formal_charge = fc;
                }

            } else if(line.compare(3,3,"END")==0){
                // End of molecule properties block
            }
        }else if(line.compare(0,3,"A  ")==0 || line.compare(0,3,"G  ")==0){
            std::getline(in, line);
        }else if(line.compare(0,3,"V  ")==0){
        }else if(line.compare(0,2,"> ")==0){
            /* If <xxx> is found on the line, use xxx as the key value */
            std::string key,val;
            size_t langle = line.find('<', 2);
            size_t rangle = line.find('>', langle+1);
            if (langle!=std::string::npos &&
                rangle!=std::string::npos) {
                key = line.substr(langle+1, rangle-langle-1);
            }
            while (std::getline(in, line)) {
                boost::trim(line);
                if (line.empty()) break;
                val=line;
            }
            add_typed_keyval(key,val,mol->ct(0));
        }else if (line.compare(0,4,"$$$$")==0){
            // End of molecule
            break;
        }else{
            /* According to the spec, we shouldnt find any other blank lines,
               but we will be forgiving just in case */
            boost::trim(line);
            if(line != "") {
                skip_to_end();
                MSYS_FAIL("Unknown line in SDF file:\n"<<line);
            }
        }
    }
    return mol;
}

SystemPtr desres::msys::ImportSdf(std::string const& path) {
    SystemPtr ct, mol = System::create();
    iterator it(path);
    while ((ct=it.next())) AppendSystem(mol,ct);
    mol->name = path;
    return mol;
}

LoadIteratorPtr desres::msys::SdfIterator(std::string const& path) {
    return LoadIteratorPtr(new iterator(path));
}

