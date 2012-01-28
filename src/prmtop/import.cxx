#include "../prmtop.hxx"
#include "../types.hxx"
#include <fstream>
#include <boost/algorithm/string.hpp>
#include <cstring>

using namespace desres::msys;

namespace {
    
    std::string parse_flag(std::string const& line) {
        std::string flag = line.substr(5);
        boost::trim(flag);
        return flag;
    }

    struct Format {
        int nperline;
        int width;
        char type;

        Format() : nperline(), width(), type() {}
        Format(std::string const& line) {
            size_t lp = line.find('(', 7);
            size_t rp = line.find(')', lp);
            if (lp==std::string::npos || rp==std::string::npos) {
                MSYS_FAIL("Expected %FORMAT(fmt), got '" << line << "'");
            }
            if (sscanf(line.c_str()+lp+1, "%d%c%d", 
                       &nperline, &type, &width)!=3) {
                MSYS_FAIL("Error parsing FORMAT '" << line << "'");
            }
        }
    };

    enum Pointers {
        Natom, Ntypes, Nbonh, Nbona, Ntheth,
        Ntheta,Nphih,  Nphia, Jparm, Nparm,
        Nnb,   Nres,   Mbona, Mtheta,Mphia,
        Numbnd,Numang, Mptra, Natyp, Nphb,
        Ifpert,Nbper,  Ngper, Ndper, Mbper,
        Mgper, Mdper,  IfBox, Nmxrs, IfCap,
        NUM_POINTERS
    };

    struct Section {
        String flag;
        Format fmt;
        String data;
    };
    typedef std::map<std::string, Section> SectionMap;

    void parse_ints(SectionMap const& map, String const& name, std::vector<int>& v) 
    {
        SectionMap::const_iterator it=map.find(name);
        if (it==map.end()) MSYS_FAIL("Missing section " << name);
        Section const& sec = it->second;
        printf("parsing %d ints for section %s\n",
               (int)v.size(), sec.flag.c_str());
        for (unsigned i=0; i<v.size(); i++) {
            std::string elem = sec.data.substr(i*sec.fmt.width, sec.fmt.width);
            if (sscanf(elem.c_str(), "%d", &v[i])!=1) {
                MSYS_FAIL("Parsing ints for " << sec.flag);
            }
        }
    }

    void parse_strs(SectionMap const& map, String const& name, 
                    std::vector<String>& v) 
    {
        SectionMap::const_iterator it=map.find(name);
        if (it==map.end()) MSYS_FAIL("Missing section " << name);
        Section const& sec = it->second;
        printf("parsing %d ints for section %s\n",
               (int)v.size(), sec.flag.c_str());
        for (unsigned i=0; i<v.size(); i++) {
            v[i] = sec.data.substr(i*sec.fmt.width, sec.fmt.width);
            boost::trim(v[i]);
        }
    }
    
}


SystemPtr desres::msys::ImportPrmTop( std::string const& path ) {

    std::string line, flag;
    std::ifstream in(path);
    if (!in) MSYS_FAIL("Could not open prmtop file at '" << path << "'");

    SystemPtr mol = System::create();
    
    /* first line is version */
    std::getline(in, line);

    /* prime the pump by looking for a %FLAG line */
    while (std::getline(in, line)) {
        if (line.size()>6 && line.substr(0,5)=="%FLAG") break;
    }
    /* slurp in the rest of the file assuming %FLAG and %FORMAT are always
       on their own line.
    */
    SectionMap section; 
    while (in) {
        std::string flag = parse_flag(line);
        Section& sec = section[flag];
        sec.flag = flag;
        while (std::getline(in, line) && line.size()<1);
        sec.fmt = Format(line);
        while (std::getline(in, line)) {
            if (line.size()<1) continue;
            if (line.substr(0,5)=="%FLAG") break;
            sec.data += line;
        }
    }

    /* build a single chain for all residues */
    Id chn = mol->addChain();

    /* build residues and atoms */
    std::vector<int> ptrs(NUM_POINTERS);
    parse_ints(section, "POINTERS", ptrs);

    std::vector<int> resptrs(ptrs[Nres]);
    parse_ints(section, "RESIDUE_POINTER", resptrs);
   
    std::vector<String> resnames(ptrs[Nres]);
    parse_strs(section, "RESIDUE_LABEL", resnames);

    std::vector<String> names(ptrs[Natom]);
    parse_strs(section, "ATOM_NAME", names);

    Id res=BadId;
    for (int i=0; i<ptrs[Natom]; i++) {
        if (i+1==resptrs[mol->residueCount()]) {
            res = mol->addResidue(chn);
            mol->residue(res).resid = mol->residueCount();
            mol->residue(res).name = resnames.at(res);
        }
        Id atm = mol->addAtom(res);
        mol->atom(atm).name = names.at(atm);
    }
    

    return mol;
}
 
