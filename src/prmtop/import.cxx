#include "../prmtop.hxx"
#include "../types.hxx"
#include "../schema.hxx"
#include "../clone.hxx"
#include "../elements.hxx"

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
        Numbnd,Numang, Nptra, Natyp, Nphb,
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
        for (unsigned i=0; i<v.size(); i++) {
            v[i] = sec.data.substr(i*sec.fmt.width, sec.fmt.width);
            boost::trim(v[i]);
        }
    }

    void parse_flts(SectionMap const& map, String const& name, std::vector<double>& v) 
    {
        SectionMap::const_iterator it=map.find(name);
        if (it==map.end()) MSYS_FAIL("Missing section " << name);
        Section const& sec = it->second;
        for (unsigned i=0; i<v.size(); i++) {
            std::string elem = sec.data.substr(i*sec.fmt.width, sec.fmt.width);
            if (sscanf(elem.c_str(), "%lf", &v[i])!=1) {
                MSYS_FAIL("Parsing dbls for " << sec.flag);
            }
        }
    }
    
}

typedef std::set<std::pair<Id,Id> > PairSet;

static void parse_nonbonded(SystemPtr mol, SectionMap const& map, int ntypes,
                            PairSet const& pairs) {

    TermTablePtr nb = AddNonbonded(mol, "vdw_12_6", "arithmetic/geometric");
    TermTablePtr pt = AddTable(mol, "pair_12_6_es");

    int ntypes2 = (ntypes * (ntypes+1))/2;
    std::vector<int> inds(ntypes*ntypes), types(mol->atomCount());
    std::vector<double> acoef(ntypes2), bcoef(ntypes2);
 
    parse_ints(map, "ATOM_TYPE_INDEX", types);
    parse_ints(map, "NONBONDED_PARM_INDEX", inds);
    parse_flts(map, "LENNARD_JONES_ACOEF", acoef);
    parse_flts(map, "LENNARD_JONES_BCOEF", bcoef);

    for (Id i=0; i<mol->atomCount(); i++) {
        int atype = types.at(i);
        int ico = inds.at((ntypes * (atype-1)+atype)-1);
        double c12 = acoef.at(ico-1);
        double c6  = bcoef.at(ico-1);
        double sig=0, eps=0;
        if (c12!=0 && c6!=0) {
            sig=pow(c12/c6, 1./6.);
            eps=c6*c6/(4*c12);
        }
        Id param = nb->params()->addParam();
        nb->params()->value(param, "sigma") = sig;
        nb->params()->value(param, "epsilon") = eps;
        IdList ids(1, i);
        nb->addTerm(ids, param);
    }
    
    for (PairSet::const_iterator it=pairs.begin(); it!=pairs.end(); ++it) {
        Id ai = it->first;
        Id aj = it->second;
        int itype = types.at(ai);
        int jtype = types.at(aj);
        int ico = inds.at((ntypes * (itype-1) + jtype)-1);
        double c12 = acoef.at(ico-1);
        double c6  = bcoef.at(ico-1);
        /* FIXME: get ES scale_factor from SCEE_SCALE_FACTOR */
        /* FIXME: get LJ scale factor from SCNB_SCALE_FACTOR */
        double lj = 1/2.0;
        double es = 1/1.2;
        double aij = lj * c12;
        double bij = lj * c6;
        double qij = es*mol->atom(ai).charge*mol->atom(aj).charge;
        Id param = pt->params()->addParam();
        pt->params()->value(param, "aij") = aij;
        pt->params()->value(param, "bij") = bij;
        pt->params()->value(param, "qij") = qij;
        IdList ids(2);
        ids[0] = ai;
        ids[1] = aj;
        pt->addTerm(ids, param);
    }
}

static void parse_stretch(SystemPtr mol, SectionMap const& map,
                          int ntypes, int nbonh, int nbona) {
    
    std::vector<double> r0(ntypes), fc(ntypes);
    std::vector<int> bonh(nbonh*3), bona(nbona*3);

    parse_flts(map, "BOND_EQUIL_VALUE", r0);
    parse_flts(map, "BOND_FORCE_CONSTANT", fc);
    parse_ints(map, "BONDS_INC_HYDROGEN", bonh);
    parse_ints(map, "BONDS_WITHOUT_HYDROGEN", bona);

    TermTablePtr tb = AddTable(mol, "stretch_harm");
    for (int i=0; i<ntypes; i++) {
        Id param = tb->params()->addParam();
        tb->params()->value(param, "fc") = fc[i];
        tb->params()->value(param, "r0") = r0[i];
    }
    IdList ids(2);
    for (int i=0; i<nbonh; i++) {
        ids[0] = bonh[3*i  ]/3;
        ids[1] = bonh[3*i+1]/3;
        tb->addTerm(ids, bonh[3*i+2]-1);
        mol->addBond(ids[0], ids[1]);
    }
    for (int i=0; i<nbona; i++) {
        ids[0] = bona[3*i  ]/3;
        ids[1] = bona[3*i+1]/3;
        tb->addTerm(ids, bona[3*i+2]-1);
        mol->addBond(ids[0], ids[1]);
    }
}
 
static void parse_angle(SystemPtr mol, SectionMap const& map,
                          int ntypes, int nbonh, int nbona) {
    
    std::vector<double> r0(ntypes), fc(ntypes);
    std::vector<int> bonh(nbonh*4), bona(nbona*4);

    parse_flts(map, "ANGLE_EQUIL_VALUE", r0);
    parse_flts(map, "ANGLE_FORCE_CONSTANT", fc);
    parse_ints(map, "ANGLES_INC_HYDROGEN", bonh);
    parse_ints(map, "ANGLES_WITHOUT_HYDROGEN", bona);

    TermTablePtr tb = AddTable(mol, "angle_harm");
    for (int i=0; i<ntypes; i++) {
        Id param = tb->params()->addParam();
        tb->params()->value(param, "fc") = fc[i];
        tb->params()->value(param, "theta0") = r0[i] * 180 / M_PI;
    }
    IdList ids(3);
    for (int i=0; i<nbonh; i++) {
        ids[0] = bonh[4*i  ]/3;
        ids[1] = bonh[4*i+1]/3;
        ids[2] = bonh[4*i+2]/3;
        tb->addTerm(ids, bonh[4*i+3]-1);
    }
    for (int i=0; i<nbona; i++) {
        ids[0] = bona[4*i  ]/3;
        ids[1] = bona[4*i+1]/3;
        ids[2] = bona[4*i+2]/3;
        tb->addTerm(ids, bona[4*i+3]-1);
    }
} 


namespace {
    struct CompareTorsion {
        bool operator() (IdList const& a, IdList const& b) const {
            for (int i=0; i<4; i++) {
                if (a[i]!=b[i]) return a[i]<b[i];
            }
            return false;
        }
    };
}
    
static void merge_torsion(ParamTablePtr params, Id id, 
                          double phase_in_radians, double fc, double period) {

    double phase = phase_in_radians * 180 / M_PI;
    /* Amber files approximate pi by 3.141594 */
    if (fabs(phase)>179.9) {
             phase = 180;
        if (phase_in_radians<0) phase *= -1;
    }
    if (phase==0) {
        /* great, nothing to do */
    } else if (phase==180) {
        /* just invert the term */
        fc *= -1;
    } else if (phase != params->value(id, "phi0").asFloat()) {
        MSYS_FAIL("multiple dihedral term contains conflicting multiplicity");
    } else {
        params->value(id, "phi0") = phase;
    }
    double oldval = params->value(id, 1+period);
    if (oldval==0) {
        params->value(id, 1+period) = fc;
    } else if (oldval != fc) {
        MSYS_FAIL("multiple dihedral term contains conflicting force constant for period " << period);
    }
    double oldsum = params->value(id, 1);
    params->value(id, 1) = oldsum + fabs(fc);
}
 
static PairSet parse_torsion(SystemPtr mol, SectionMap const& map,
                          int ntypes, int nbonh, int nbona) {
    
    std::vector<double> phase(ntypes), fc(ntypes), period(ntypes);
    std::vector<int> bonh(nbonh*5), bona(nbona*5);

    parse_flts(map, "DIHEDRAL_PHASE", phase);
    parse_flts(map, "DIHEDRAL_FORCE_CONSTANT", fc);
    parse_flts(map, "DIHEDRAL_PERIODICITY", period);
    parse_ints(map, "DIHEDRALS_INC_HYDROGEN", bonh);
    parse_ints(map, "DIHEDRALS_WITHOUT_HYDROGEN", bona);

    bonh.insert(bonh.end(), bona.begin(), bona.end());
    bona.clear();

    /* hash the torsion terms, converting negative indices to positive as
       needed and tracking which terms should generate 1-4 pair terms.  */

    PairSet pairs;
    typedef std::map<IdList, Id, CompareTorsion> TorsionHash;
    TorsionHash hash;
    IdList ids(4);
    TermTablePtr tb = AddTable(mol, "dihedral_trig");
    
    for (int i=0; i<nbonh+nbona; i++) {
        int ai = bonh[5*i  ]/3;
        int aj = bonh[5*i+1]/3;
        int ak = bonh[5*i+2]/3;
        int al = bonh[5*i+3]/3;
        int ind = bonh[5*i+4]-1;
        bool needs_pair = false;
        if (ak<0) {
            ak = -ak;
        } else {
            needs_pair = true;
        }
        if (al<0) { /* it's an improper, though we treat it the same */
            al = -al;
        }
        ids[0]=ai;
        ids[1]=aj;
        ids[2]=ak;
        ids[3]=al;
        std::pair<TorsionHash::iterator, bool> r;
        r = hash.insert(std::make_pair(ids, BadId ));
        Id param = BadId;
        if (r.second) {
            /* new term, so create a Term and Param entry */
            param = tb->params()->addParam();
            tb->addTerm(ids, param);
            r.first->second = param;
        } else {
            param = r.first->second;
        }
        if (needs_pair) {
            Id pi=ai, pj=al;
            if (pi>pj) std::swap(pi,pj);
            pairs.insert(std::make_pair(pi,pj));
        }
        merge_torsion(tb->params(), param, phase.at(ind), fc.at(ind), period.at(ind));
    }
    return pairs;
} 

static void parse_exclusions(SystemPtr mol, SectionMap const& map, int n) {
    if (n==0) return;
    TermTablePtr tb=AddTable(mol, "exclusion");
    std::vector<int> nexcl(mol->atomCount()), excl(n);
    parse_ints(map, "NUMBER_EXCLUDED_ATOMS", nexcl);
    parse_ints(map, "EXCLUDED_ATOMS_LIST", excl);
    unsigned j=0;
    IdList ids(2);
    for (Id ai=0; ai<mol->atomCount(); ai++) {
        ids[0]=ai;
        for (int i=0; i<nexcl[ai]; i++, j++) {
            Id aj = excl[j];
            if (aj==0) continue;
            ids[1]=aj-1;
            tb->addTerm(ids, BadId);
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

    std::vector<Float> charges(ptrs[Natom]);
    parse_flts(section, "CHARGE", charges);

    std::vector<Float> masses(ptrs[Natom]);
    parse_flts(section, "MASS", masses);

    Id res=BadId;
    resptrs.push_back(ptrs[Natom]+1); /* simplify residue start logic */
    for (int i=0; i<ptrs[Natom]; i++) {
        if (i+1==resptrs[mol->residueCount()]) {
            res = mol->addResidue(chn);
            mol->residue(res).resid = mol->residueCount();
            mol->residue(res).name = resnames.at(res);
        }
        Id atm = mol->addAtom(res);
        mol->atom(atm).name = names.at(atm);
        mol->atom(atm).charge = charges.at(atm) / 18.2223; /* magic scale */
        mol->atom(atm).mass = masses.at(atm);
        mol->atom(atm).atomic_number = GuessAtomicNumber(masses.at(atm));
    }

    parse_stretch(mol, section, ptrs[Numbnd], ptrs[Nbonh], ptrs[Nbona]);
    parse_angle(mol, section, ptrs[Numang], ptrs[Ntheth], ptrs[Ntheta]);
    PairSet pairs = parse_torsion(mol, section, 
                                  ptrs[Nptra], ptrs[Nphih], ptrs[Nphia]);
    parse_nonbonded(mol, section, ptrs[Ntypes], pairs);
    parse_exclusions(mol, section, ptrs[Nnb]);

    mol->updateFragids();
    mol->coalesceTables();
    return Clone(mol, mol->atoms());
}
 
