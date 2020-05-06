#include "token.hxx"
#include "../elements.hxx"
#include "../smarts.hxx"
#include <unordered_map>
#include <unordered_set>
//#include <regex>
#include <boost/regex.hpp>

using namespace desres::msys;
using namespace desres::msys::atomsel;

Selection desres::msys::atomsel::full_selection(System* sys) {
    Selection s(sys->maxAtomId());
    s.fill();
    if (sys->maxAtomId() != sys->atomCount()) {
        for (Id i=0, n=sys->maxAtomId(); i<n; i++) {
            if (!sys->hasAtom(i)) s[i]=0;
        }
    }
    return s;
}

static int get_index(Query* q, Id i) {
    return i;
}
static int get_frag(Query* q, Id i) {
    auto mol = q->mol;
    return mol->atomFAST(i).fragid;
}
static int get_numbonds(Query* q, Id i) {
    auto mol = q->mol;
    return mol->bondCountForAtom(i);
}
static int get_degree(Query* q, Id i) {
    auto mol = q->mol;
    if (mol->atomFAST(i).atomic_number == 0) return 0;
    int degree = 0;
    for (auto bid : mol->bondsForAtom(i)) {
        bond_t const& bnd = mol->bondFAST(bid);
        if (mol->atomFAST(bnd.other(i)).atomic_number > 0) {
            ++degree;
        }
    }
    return degree;
}

static int get_ctnumber(Query* q, Id i) {
    auto mol = q->mol;
    return mol->chainFAST(mol->residueFAST(mol->atomFAST(i).residue).chain).ct+1;
}
static int get_ct(Query* q, Id i) {
    auto mol = q->mol;
    return mol->chainFAST(mol->residueFAST(mol->atomFAST(i).residue).chain).ct;
}
static int get_resid(Query* q, Id i) {
    auto mol = q->mol;
    return mol->residueFAST(mol->atomFAST(i).residue).resid;
}
static int get_residue(Query* q, Id i) {
    auto mol = q->mol;
    return mol->atomFAST(i).residue;
}
static int get_anum(Query* q, Id i) {
    auto mol = q->mol;
    return mol->atomFAST(i).atomic_number;
}
static double get_mass(Query* q, Id i) {
    auto mol = q->mol;
    return mol->atomFAST(i).mass;
}
static double get_charge(Query* q, Id i) {
    auto mol = q->mol;
    return mol->atomFAST(i).charge;
}
static double get_x(Query* q, Id i) {
    return q->pos ? (double)q->pos[3*i] : q->mol->atomFAST(i).x;
}
static double get_y(Query* q, Id i) {
    return q->pos ? (double)q->pos[3*i+1] : q->mol->atomFAST(i).y;
}
static double get_z(Query* q, Id i) {
    return q->pos ? (double)q->pos[3*i+2] : q->mol->atomFAST(i).z;
}
static double get_vx(Query* q, Id i) {
    return q->mol->atomFAST(i).vx;
}
static double get_vy(Query* q, Id i) {
    return q->mol->atomFAST(i).vy;
}
static double get_vz(Query* q, Id i) {
    return q->mol->atomFAST(i).vz;
}
static const char* get_chain(Query* q, Id i) {
    auto mol = q->mol;
    return mol->chainFAST(mol->residueFAST(mol->atomFAST(i).residue).chain).name.data();
}
static const char* get_segid(Query* q, Id i) {
    auto mol = q->mol;
    return mol->chainFAST(mol->residueFAST(mol->atomFAST(i).residue).chain).segid.c_str();
}
static const char* get_name(Query* q, Id i) {
    auto mol = q->mol;
    return mol->atomFAST(i).name.c_str();
}
static const char* get_resname(Query* q, Id i) {
    auto mol = q->mol;
    return mol->residueFAST(mol->atomFAST(i).residue).name.c_str();
}
static const char* get_insertion(Query* q, Id i) {
    auto mol = q->mol;
    return mol->residueFAST(mol->atomFAST(i).residue).insertion.c_str();
}
static const char* get_element(Query* q, Id i) {
    auto mol = q->mol;
    return AbbreviationForElement(mol->atomFAST(i).atomic_number);
}
static int get_iprop(Query* q, Id i) {
    auto mol = q->mol;
    return mol->atomPropValue(i, q->id).asInt();
}
static double get_fprop(Query* q, Id i) {
    auto mol = q->mol;
    return mol->atomPropValue(i, q->id).asFloat();
}
static const char*get_sprop(Query* q, Id i) {
    auto mol = q->mol;
    return mol->atomPropValue(i, q->id).c_str();
}

typedef int (*iget)(Query*,Id);
typedef double (*fget)(Query*,Id);
typedef const char* (*sget)(Query*,Id);

enum {
    StrFuncType = 4
};
struct Getter {
    const int type;
    union {
        iget i;
        fget f;
        sget s;
    };
    Getter(iget i) : type(IntType), i(i) {}
    Getter(fget f) : type(FloatType), f(f) {}
    Getter(sget s) : type(StringType), s(s) {}
};

static const std::unordered_map<std::string, Getter> map = {
    {"atomicnumber", get_anum},
    {"mass", get_mass},
    {"charge", get_charge},
    {"chain", get_chain},
    {"ctnumber", get_ctnumber},
    {"ct", get_ct},
    {"element", get_element},
    {"fragment", get_frag},
    {"fragid", get_frag},
    {"index", get_index},
    {"name", get_name},
    {"numbonds", get_numbonds},
    {"degree", get_degree},
    {"resid", get_resid},
    {"residue", get_residue},
    {"resname", get_resname},
    {"insertion", get_insertion},
    {"segid", get_segid},
    {"segname", get_segid},
    {"x", get_x},
    {"y", get_y},
    {"z", get_z},
    {"vx", get_vx},
    {"vy", get_vy},
    {"vz", get_vz},
};

template <typename T>
static bool find_in_ranges(std::vector<std::pair<T,T>> const& rng, T i) {
    for (auto const& r : rng) {
        if (i>=r.first && i<=r.second) return true;
    }
    return false;
}

static void eval_smarts(System* mol, std::vector<std::string> const& pats,
        Selection& s) {
    AnnotatedSystem a(mol->shared_from_this());
    /* get fragids of selected atoms */
    std::unordered_set<Id> fragids;
    for (Id i=0, n=s.size(); i<n; i++) {
        if (s[i]) fragids.insert(mol->atomFAST(i).fragid);
    }
    /* use as starts atoms in same fragment as selected atoms */
    IdList starts;
    for (Id i=0, n=s.size(); i<n; i++) {
        if (fragids.count(mol->atomFAST(i).fragid)) starts.push_back(i);
    }

    /* find atoms matched by smarts pattern */
    Selection sel(s.size());
    for (auto&& pat : pats) {
        SmartsPattern p(pat);
        for (auto&& ids : p.findMatches(a, starts)) {
            //printf("hits: %u\n", (Id)ids.size());
            for (auto id : ids) {
                sel[id]=1;
            }
        }
    }
    s.intersect(sel);
}

static void eval_paramtype(System* mol, std::vector<std::string> const& pats,
        Selection& s) {
    if (pats.empty()) {
        MSYS_FAIL("paramtype selection requires table name as first argument");
    }
    auto table = mol->table(pats[0]);
    if (!table) {
        MSYS_FAIL("paramtype selection references nonexistent table '" << pats[0] << "'");
    }
    auto params = table->params();
    Id col = params->propIndex("type");
    if (bad(col)) {
        MSYS_FAIL("paramtype selection references table '" << pats[0] << "' with no 'type' column");
    }
    std::unordered_set<Id> ids;
    for (Id i=1; i<pats.size(); i++) {
        auto tmp = params->findString(col, pats[i]);
        ids.insert(tmp.begin(), tmp.end());
    }
    const Id natoms = table->atomCount();
    Selection sel(mol->atomCount());
    for (auto i=table->begin(), e=table->end(); i!=e; ++i) {
        if (ids.count(i->param())) {
            for (Id j=0; j<natoms; j++) {
                sel[i->atom(j)] = 1;
            }
        }
    }
    s.intersect(sel);
}

static std::unordered_map<std::string, char> single_letter_residue_map({
        { "GLY", 'G'},
        { "ALA", 'A'},
	{ "VAL", 'V'},
	{ "PHE", 'F'},
	{ "PRO", 'P'},
	{ "MET", 'M'},
	{ "ILE", 'I'},
	{ "LEU", 'L'},
	{ "ASP", 'D'},
	{ "GLU", 'E'},
	{ "LYS", 'K'},
	{ "ARG", 'R'},
	{ "SER", 'S'},
	{ "THR", 'T'},
	{ "TYR", 'Y'},
	{ "HIS", 'H'},
	{ "CYS", 'C'},
	{ "ASN", 'N'},
	{ "GLN", 'Q'},
	{ "TRP", 'W'},
        // HIS synonyms
        { "HSE", 'H'},
        { "HSD", 'H'},
        { "HSP", 'H'},
        { "HID", 'H'},
        { "HIP", 'H'},
        // CYS synonyms
        { "CYX", 'C'},
        // LYS synonyms
        { "LYP", 'K'},
        // nucleics
        { "ADE", 'A'},
        { "A", 'A'},
        { "THY", 'T' },
        { "T", 'T' },
        { "CYT", 'C' },
        { "C", 'C' },
        { "GUA", 'G' },
        { "G", 'G' },
});

static void eval_sequence(System* mol, std::vector<std::string> const& pats, Selection& s) {

    std::vector<boost::regex> regs;
    for (auto const& pat : pats) {
        regs.emplace_back(pat);
    }
    // convert each chain to a list of single character residue names.
    // match the regex against the chain.
    // for every hit, turn on the corresponding atoms in the residues
    Selection sel(mol->maxAtomId());
    for (auto chain_id : mol->chains()) {
        std::string codes;
        for (auto residue_id : mol->residuesForChain(chain_id)) {
            auto const& res = mol->residueFAST(residue_id);
            if (res.type == ResidueWater) {
                codes.push_back('X');
            } else {
                auto p = single_letter_residue_map.find(mol->residueFAST(residue_id).name);
                char c = (p == single_letter_residue_map.end()) ? 'X' : p->second;
                codes.push_back(c);
            }
        }
        for (auto const& reg : regs) {
            auto match = boost::sregex_iterator(codes.begin(), codes.end(), reg);
            auto end = boost::sregex_iterator();
            for (; match!=end; ++match) {
                for (auto i=match->position(), j=match->position() + match->length(); i!=j; ++i) {
                    auto residue_id = mol->residuesForChain(chain_id).at(i);
                    for (auto atmid : mol->atomsForResidue(residue_id)) {
                        sel[atmid] = 1;
                    }
                }
            }
        }
    }
    s.intersect(sel);
}


typedef void (*sfunc)(System*, std::vector<std::string>const&, Selection&);
static const std::unordered_map<std::string,sfunc> strfuncs = {
    {"smarts", eval_smarts},
    {"paramtype", eval_paramtype},
    {"sequence", eval_sequence},
};

bool desres::msys::atomsel::is_keyword(std::string const& name, System* mol) {
    return map.find(name) != map.end()
        || strfuncs.find(name) != strfuncs.end()
        || mol->atomPropIndex(name)!=BadId
        ;
}

static Getter lookup(std::string const& name, Query* q) {
    auto iter = map.find(name);
    if (iter!=map.end()) {
        return iter->second;
    }
    Id col = q->mol->atomPropIndex(name);
    if (bad(col)) {
        MSYS_FAIL("Unknown atom selection keyword '" << name << "'");
    }
    q->id = col;
    switch (q->mol->atomPropType(col)) {
        case IntType:   return get_iprop;
        case FloatType: return get_fprop;
        default:;
        case StringType: return get_sprop;
    };
}


void KeyExpr::eval(Selection const& s, std::vector<double>& v) {
    auto g = lookup(name, q);
    switch (g.type) {
        case IntType: 
            for (Id i=0, n=s.size(); i<n; i++) if (s[i]) v[i] = g.i(q, i);
            break;
        case FloatType:
            for (Id i=0, n=s.size(); i<n; i++) if (s[i]) v[i] = g.f(q, i);
            break;
        default:;
        case StringType:
            MSYS_FAIL("selection keyword '" << name << "' with string type cannot be used in a numeric expression");
    }
}

void KeyPredicate::eval(Selection& s) {
    auto f = strfuncs.find(name);
    if (f!=strfuncs.end()) {
        f->second(q->mol, va->sval, s);
        return;
    }

    auto g = lookup(name, q);
    std::vector<boost::regex> regs;
    switch (g.type) {
    case IntType:
        {
            if (!(va->fval.empty() && va->sval.empty() && va->frng.empty() && va->regex.empty())) {
                MSYS_FAIL("Selection keyword '" << name << "' expects integer values, but got float, string or regex values.");
            }
            sort_unique(va->ival);
            for (Id i=0, n=s.size(); i<n; i++) {
                if (s[i]) {
                    int ival = g.i(q,i);
                    if (std::binary_search(va->ival.begin(), va->ival.end(),ival)) continue;
                    if (find_in_ranges(va->irng, ival)) continue;
                    s[i] = 0;
                }
            }
        }
        break;
    case FloatType:
            if (!(va->sval.empty() && va->regex.empty())) {
                MSYS_FAIL("Selection keyword '" << name << "' expects integer values, but got string or regex values.");
            }
            sort_unique(va->ival); 
            sort_unique(va->fval); 
            for (Id i=0, n=s.size(); i<n; i++) {
                if (s[i]) {
                    double val = g.f(q,i);
                    if (std::binary_search(va->ival.begin(), va->ival.end(),(int)val)) continue;
                    if (std::binary_search(va->fval.begin(), va->fval.end(),val)) continue;
                    if (find_in_ranges(va->irng, (int)val)) continue;
                    if (find_in_ranges(va->frng, val)) continue;
                    s[i] = 0;
                }
            }
            break;
    case StringType:
            if (!(va->ival.empty() && va->fval.empty() && va->irng.empty() && va->frng.empty())) {
                MSYS_FAIL("Selection keyword '" << name << "' expects string values, but got float or integer values.");
            }
            sort_unique(va->sval);
            for (auto const& r : va->regex) {
                regs.emplace_back(r);
            }
            for (Id i=0, n=s.size(); i<n; i++) {
                if (s[i]) {
                // TODO: compare the values as const char* rather than 
                // constructing std::string intermediates
                    std::string val(g.s(q,i));
                    if (std::binary_search(va->sval.begin(), va->sval.end(), std::string(val))) continue;
                    boost::smatch match;
                    bool found = false;
                    for (auto const& re : regs) {
                        if (boost::regex_match(val, match, re)) {
                            found = true;
                            break;
                        }
                    }
                    if (found) continue;
                    s[i]=0;
                }
            }
            break;
    case StrFuncType:
            /* FIXME: can never get here because of logic at the top of
             * this function */
            if (!(va->ival.empty() && va->fval.empty() && va->irng.empty() && va->frng.empty() && va->regex.empty())) {
                MSYS_FAIL("Selection keyword '" << name << "' expects only string arguments, but got float, integer or regex.");
            }


    default:;
    }
}


void SamePredicate::eval(Selection& s) {
    Selection x = full_selection(q->mol);
    sub->eval(x);

    auto g = lookup(name, q);
    switch (g.type) {
        case IntType: 
        {
            std::unordered_set<int> h;
            for (Id i=0, n=x.size(); i<n; i++) {
                if (x[i]) h.insert(g.i(q,i));
            }
            for (Id i=0, n=s.size(); i<n; i++) {
                if (s[i]) s[i] = h.count(g.i(q,i));
            }
        }
        break;
        case FloatType:
        {
            std::unordered_set<double> h;
            for (Id i=0, n=x.size(); i<n; i++) {
                if (x[i]) h.insert(g.f(q,i));
            }
            for (Id i=0, n=s.size(); i<n; i++) {
                if (s[i]) s[i] = h.count(g.f(q,i));
            }
        }
        break;
        case StringType:
        {
            std::unordered_set<std::string> h;
            for (Id i=0, n=x.size(); i<n; i++) {
                if (x[i]) h.insert(g.s(q,i));
            }
            for (Id i=0, n=s.size(); i<n; i++) {
                if (s[i]) s[i] = h.count(g.s(q,i));
            }
        }
        break;
    }
}
