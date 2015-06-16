#include "../sdf.hxx"
#include "elements.hxx"
#include "../append.hxx"
#include <boost/foreach.hpp>
#include <boost/algorithm/string.hpp> /* for boost::trim */
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/math/special_functions/fpclassify.hpp>

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
            if (boost::math::isfinite(v)) {
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
        explicit iterator(std::istream& in) {
            init(in);
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

    class scanner : public MoleculeIterator {
        FILE* fp = 0;
        char buf[82];
        bool getline() {
            return fgets(buf, sizeof(buf), fp)!=NULL;
        }
        bool eof() {
            return feof(fp);
        }
        std::string skip_to_end() {
            std::string current(buf);
            while (getline()) {
                if (!strncmp(buf, "$$$$", 4)) break;
            }
            return current;
        }
        // parse text of the form __z, _yz, or xyz, where underscore
        // denotes leading space and x,y,z are digits.  Return -1 on
        // error.
        static short parse_count(const char* s) {
            int n=0;
            unsigned short result=0;
            unsigned short sign = 1;
            if (isspace(*s)) { ++n; ++s; }
            if (isspace(*s)) { ++n; ++s; }
            if (*s == '-') {
                sign = -sign;
                ++s;
                ++n;
            }
            for (; n<3; n++, s++) {
                if (!isdigit(*s)) return -1;
                unsigned short d = *s - '0';
                switch (n) {
                    case 0: result += 100*d; break;
                    case 1: result +=  10*d; break;
                    case 2: result +=     d; break;
                };
            }
            return result;
        }

        // parse xxxxx.yyyy as float
        static float parse_coord(const char* s) {
            float scale = 10000.0f;
            unsigned result = 0;
            static const float fail = std::numeric_limits<float>::max();
            int n=0;
            for (; n<4; n++, s++) {
                if (!isspace(*s)) break;
            }
            if (*s == '-') {
                ++s;
                ++n;
                scale = -scale;
            }
            for (; n<5; n++, s++) {
                if (!isdigit(*s)) return fail;
                unsigned d = *s - '0';
                switch (n) {
                    case 0: result += 100000000*d; break;
                    case 1: result +=  10000000*d; break;
                    case 2: result +=   1000000*d; break;
                    case 3: result +=    100000*d; break;
                    case 4: result +=     10000*d; break;
                }
            }
            if (*s != '.') return fail;
            result += 1000*(s[1]-'0');
            result +=  100*(s[2]-'0');
            result +=   10*(s[3]-'0');
            result +=      (s[4]-'0');
            return float(result)/scale;
        }

        static char parse_element(const char* s) {
            if (isspace(*s)) ++s;
            if (isspace(*s)) {
                // one-character element
                ++s;
                switch (*s) {
                    case 'C': return 6;
                    case 'H': return 1;
                    case 'N': return 7;
                    case 'O': return 8;
                    case 'F': return 9;
                    case 'P': return 15;
                    case 'S': return 16;
                    case 'K': return 19;
                    case 'B': return 5;
                    case 'V': return 23;
                    case 'Y': return 39;
                    case 'I': return 53;
                    case 'W': return 74;
                    case 'U': return 92;
                    default:  return 0;
                };
            }
            char buf[3] = {s[0],s[1],'\0'};
            return ElementForAbbreviationSlow(buf);
        }

    public:
        ~scanner() { 
            if (fp) fclose(fp); 
        }
        explicit scanner(std::string path) {
            fp = fopen(path.c_str(), "r");
            if (!fp) {
                MSYS_FAIL("Error opening " << path << " for reading: " 
                        << strerror(errno));
            }
        }
        MoleculePtr next() {
            MoleculePtr ptr;

            // three lines for header block, then one for counts
            getline();
            std::string name = buf;
            name.pop_back();    // remove trailing newline
            getline();
            getline();
            getline();
            if (feof(fp)) return ptr;
            auto natoms = parse_count(buf);
            auto nbonds = parse_count(buf+3);
            if (bad(natoms) || bad(nbonds)) {
                MSYS_FAIL("Bad counts line: " << skip_to_end());
            }

            ptr.reset(new Molecule(natoms, nbonds));
            ptr->name() = name;

            // atoms
            for (unsigned short i=0; i<natoms; i++) {
                if (!getline()) {
                    skip_to_end();
                    MSYS_FAIL("Missing expected atom record");
                }
                auto sz = strlen(buf);
                if (sz<34) {
                    MSYS_FAIL("Malformed atom line: " << skip_to_end());
                }
                ptr->atom(i).x = parse_coord(buf   );
                ptr->atom(i).y = parse_coord(buf+10);
                ptr->atom(i).z = parse_coord(buf+20);
                ptr->atom(i).atomic_number = parse_element(buf+31);

                if (sz>=39) {
                    auto q = parse_count(buf+36);
                    ptr->atom(i).formal_charge = q==0 ? 0 : 4-q;
                }
                if (sz>=42) {
                    auto s = parse_count(buf+39);
                    ptr->atom(i).stereo_parity = s;
                }
            }

            // bonds
            for (unsigned short i=0; i<nbonds; i++) {
                if (!getline()) {
                    skip_to_end();
                    MSYS_FAIL("Missing expected bond record");
                }
                auto sz = strlen(buf);
                if (sz<9) {
                    MSYS_FAIL("Malformed bond line: " << skip_to_end());
                }
                ptr->bond(i).i =     parse_count(buf  )-1;
                ptr->bond(i).j =     parse_count(buf+3)-1;
                ptr->bond(i).order = parse_count(buf+6);
                ptr->bond(i).stereo = parse_count(buf+9);
            }

            // M entries
            while (getline()) {
                if (!strncmp(buf, "M  ", 3)) {
                    if (!strncmp(buf+3, "END", 3)) {
                        break;
                    } else if (!strncmp(buf+3, "CHG", 3)) {
                        for (int i=0; i<natoms; i++) {
                            ptr->atom(i).formal_charge = 0;
                        }
                        int n = parse_count(buf+6);
                        for (int i=0; i<n; i++) {
                            short aid = parse_count(buf+10+8*i);
                            short chg = parse_count(buf+14+8*i);
                            if (aid<1) {
                                MSYS_FAIL("Malformed CHG line: " << skip_to_end());
                            }
                            ptr->atom(aid-1).formal_charge = chg;
                        }
                    }
                } else if (!strncmp(buf, "A  ", 3)) {
                    getline();
                } else if (!strncmp(buf, "G  ", 3)) {
                    getline();
                } else if (!strncmp(buf, "V  ", 3)) {
                    // ignore
                } else {
                    MSYS_FAIL("Malformed properties line: " << skip_to_end());
                }
            }

            // data fields.  
            bool needline = true;
            for (;;) {
                if (needline && !getline()) {
                    MSYS_FAIL("Unexpected end of file");
                }
                needline = true;
                if (!strncmp(buf, "$$$$", 4)) {
                    break;
                } else if (!strncmp(buf, "> ", 2)) {
                    const char* langle = strchr(buf+2,'<');
                    const char* rangle = strchr(buf+3,'>');
                    std::string key;
                    if (langle && rangle) {
                        key = std::string(langle+1,rangle);
                    }
                    std::string val;
                    // support multiline values containing newlines
                    for (;;) {
                        if (!getline()) {
                            MSYS_FAIL("Unexpected end of file");
                        }
                        if (!strncmp(buf, "$$$$", 4) ||
                            !strncmp(buf, "> ", 2)) {
                            if (*key.c_str()) {
                                auto sz = val.size();
                                if (sz==1) val.clear(); /* just a newline */
                                else if (sz>1) val.resize(sz-2);
                                ptr->data().insert(
                                        Molecule::Data::value_type(key,val));
                            }
                            needline = false;
                            break;
                        }
                        val += buf;
                    }
                } else {
                    MSYS_FAIL("Malformed data field line: " << skip_to_end());
                }
            }
            return ptr;
        }
    };
}

MoleculeIteratorPtr desres::msys::ScanSdf(std::string path) {
    return std::unique_ptr<MoleculeIterator>(new scanner(path));
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
    //int nbonds = stringToInt(line.substr(3,3));
    /* HACK: Support lines like that which were produced by older msys:
9301001  0  0  1  0            999 V2000
    */
    std::string bondstr = line.substr(3,4);
    boost::trim(bondstr);
    int nbonds = stringToInt(bondstr);

    //std::vector<int> valenceList;

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
                MSYS_FAIL("Skipping entry with doublet radical charge specification:\n" << line);
            }
            /* Just set formal charge, not 'charge', for consistency with
             * the exporter, which only looks at formal charge.  */
            atom.formal_charge = q==0 ? 0 : 4-q;
            
            // we don't currently use the valence info, so let's not waste
            // our time parsing it.
#if 0
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
#endif
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
    bool needline = true;
    for(;;) {
        if (needline && !std::getline(in, line)) break;
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
                if(line.compare(0,2,"> ")==0 ||
                   line.compare(0,4,"$$$$")==0) break;
                if (!val.empty()) val += '\n';
                val += line;
            }
            boost::trim(val);
            add_typed_keyval(key,val,mol->ct(0));
            needline = false;
            continue;
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
        needline = true;
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

SystemPtr desres::msys::ImportSdfFromStream(std::istream& in) {
    SystemPtr ct, mol = System::create();
    iterator it(in);
    while ((ct=it.next())) AppendSystem(mol,ct);
    return mol;
}


LoadIteratorPtr desres::msys::SdfIterator(std::string const& path) {
    return LoadIteratorPtr(new iterator(path));
}

