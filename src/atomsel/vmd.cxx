#include "vmd.hxx"
#include "msys_keyword.hxx"
#include "vmd_keyword.hxx"
#include "keyword_predicate.hxx"
#include "within_predicate.hxx"

#include "../smarts.hxx"

using namespace desres::msys;
using namespace desres::msys::atomsel;

int vmd_parse(VMD* ctxt);

namespace {
  struct atomsel_macro {
      const char * name;
      const char * text;
  };

  atomsel_macro builtin_macros[] = {
      {"at","resname ADE A THY T"},
      {"acidic","resname ASP GLU"},
      {"cyclic","resname HIS PHE PRO TRP TYR"},
      {"acyclic","protein and not cyclic"},
      {"aliphatic","resname ALA GLY ILE LEU VAL"},
      {"alpha","protein and name CA"},
      {"amino","protein"},
      {"aromatic","resname HIS PHE TRP TYR"},
      {"basic","resname ARG HIS LYS HSP"},
      {"bonded","numbonds > 0"},
      {"buried"," resname ALA LEU VAL ILE PHE CYS MET TRP"},
      {"cg","resname CYT C GUA G"},
      {"charged","basic or acidic"},
      {"hetero","not (protein or nucleic)"},

    // APD's hydrophobic residue list, from Branden and Tooze (pp6-7).
      {"hydrophobic","resname ALA LEU VAL ILE PRO PHE MET TRP"},
    
      {"small","resname ALA GLY SER"},
      {"medium","resname VAL THR ASP ASN PRO CYS ASX PCA HYP"},
      {"large","protein and not (small or medium)"},
      {"neutral","resname VAL PHE GLN TYR HIS CYS MET TRP ASX GLX PCA HYP"},
      {"polar","protein and not hydrophobic"},
      {"purine","resname ADE A GUA G"},
      {"pyrimidine","resname CYT C THY T URA U"},
      {"surface","protein and not buried"},
      {"lipid","resname DLPE DMPC DPPC GPC LPPC PALM PC PGCL POPC POPE POPS"},
      {"lipids","lipid"},
      {"ion","resname AL BA CA Ca CAL CD CES CLA CL Cl CO CS CU Cu CU1 CUA HG IN IOD K MG MN3 MO3 MO4 MO5 MO6 NA Na NAW OC7 PB POT PT RB SOD TB TL WO4 YB ZN ZN1 ZN2"},
      {"ions","ion"},
      {"sugar","resname AGLC"},
      {"solvent","not (protein or sugar or nucleic or lipid)"},
      /* for carbon, nitrogen, oxygen, and sulfur (and hydrogen), VMD
       * uses a name-based regex; e.g. 'name "N.*"' for nitrogen.  This
       * is just silly, and gets things like Na and Cl wrong.  We refuse
       * to reproduce this buggy behavior and intead look to the atomic
       * number. */
      {"carbon","atomicnumber 6"},
      {"nitrogen","atomicnumber 7"},
      {"oxygen","atomicnumber 8"},
      {"sulfur","atomicnumber 16"},
      {"noh","not hydrogen"},
      {"heme","resname HEM HEME"}
  };
}

extern int vmd_debug;

PredicatePtr VMD::parse(std::string const& sel) {
    //vmd_debug = 1;
    txt = sel.data();
    result = NULL;
    init_scanner();
    vmd_parse(this);
    destroy_scanner();
    return result ? result->shared_from_this() : PredicatePtr();
}

Predicate* VMD::add(PredicatePtr p) {
    if (p) {
        predicates.push_back(p);
        return p.get();
    }
    return NULL;
}
Expression* VMD::add(ExpressionPtr e) {
    if (e) {
        expressions.push_back(e);
        return e.get();
    }
    return NULL;
}

Keyword* VMD::find_key(const char* id) {
    KeywordPtr key = 
      !strcmp(id,"atomicnumber") ? keyword_anum(sys) :
      !strcmp(id,"chain") ? keyword_chain(sys) :
      !strcmp(id,"segid") ? keyword_segid(sys) :
      !strcmp(id,"segname") ? keyword_segid(sys) :
      !strcmp(id,"ctnumber") ? keyword_ctnumber(sys) :
      !strcmp(id,"charge") ? keyword_charge(sys) :
      !strcmp(id,"element") ? keyword_element(sys) :
      !strcmp(id,"fragment") ? keyword_fragment(sys) :
      !strcmp(id,"index") ? keyword_index(sys) :
      !strcmp(id,"mass") ? keyword_mass(sys) :
      !strcmp(id,"name") ? keyword_name(sys) :
      !strcmp(id,"numbonds") ? keyword_numbonds(sys) :
      !strcmp(id,"resid") ? keyword_resid(sys) :
      !strcmp(id,"residue") ? keyword_residue(sys) :
      !strcmp(id,"resname") ? keyword_resname(sys) :
      !strcmp(id,"insertion") ? keyword_insertion(sys) :
      !strcmp(id,"fragid") ? keyword_fragid(sys) :
      !strcmp(id,"x") ? keyword_x(sys) :
      !strcmp(id,"y") ? keyword_y(sys) :
      !strcmp(id,"z") ? keyword_z(sys) :
      !strcmp(id,"vx") ? keyword_x(sys) :
      !strcmp(id,"vy") ? keyword_y(sys) :
      !strcmp(id,"vz") ? keyword_z(sys) :
      KeywordPtr();
    if (!key) key = keyword_atomprop(sys, id);
    if (key) {
      keywords.push_back(key);
    }
    return key.get();
}

Keyword* VMD::find_single(const char* id) {
    KeywordPtr key = 
        !strcmp(id,"water") ? keyword_water(sys) :
        !strcmp(id,"hydrogen") ? keyword_hydrogen(sys) :
        !strcmp(id,"backbone") ? keyword_backbone(sys) :
        !strcmp(id,"sidechain") ? keyword_sidechain(sys) :
        !strcmp(id,"protein") ? keyword_protein(sys) :
        //!strcmp(id,"alchemical") ? keyword_alchemical(sys) :
        !strcmp(id,"nucleic") ? keyword_nucleic(sys) :
        !strcmp(id,"all") ? keyword_all() :
        !strcmp(id,"none") ? keyword_none() :
        KeywordPtr();
    if (key) {
        keywords.push_back(key);
    }
    return key.get();
}

namespace {
    class SmartsPredicate : public StringPredicate {
        AnnotatedSystem sys;
        std::vector<SmartsPattern> pats;
    public:
        SmartsPredicate(SystemPtr sys) : sys(sys) {}
        void add(std::string const& s) {
            pats.push_back(SmartsPattern(s));
        }
        void eval(Selection& s) {
            Selection sub(s.size());
            IdList starts = sys.atoms();
            BOOST_FOREACH(SmartsPattern const& pat, pats) {
                BOOST_FOREACH(IdList ids, pat.findMatches(sys, starts)) {
                    BOOST_FOREACH(Id id, ids) {
                        sub[id]=1;
                    }
                }
            }
            s.intersect(sub);
        }
        void dump(std::ostream& out) const {}
    };

    class ParamTypePredicate : public StringPredicate {
        SystemPtr mol;
        std::string tablename;
        std::vector<std::string> types;
    public:
        ParamTypePredicate(SystemPtr mol) : mol(mol) {}
        void add(std::string const& s) {
            if (tablename.empty()) {
                tablename = s;
            } else {
                types.push_back(s);
            }
        }
        void eval(Selection& s) {
            if (tablename.empty()) {
                MSYS_FAIL("Missing table name");
            }
            TermTablePtr table = mol->table(tablename);
            if (!table) {
                MSYS_FAIL("No such table '" << tablename << "'");
            }
            ParamTablePtr params = table->params();
            Id index = params->propIndex("type");
            if (bad(index)) {
                MSYS_FAIL("Table '" << tablename << "' has no 'type' column");
            }
            /* get params matching selection criteria */
            IdList ids;
            for (auto const& type : types) {
                IdList tmp = params->findString(index, type);
                ids.insert(ids.end(), tmp.begin(), tmp.end());
            }
            sort_unique(ids);
            /* get atoms in terms with these params */
            const Id natoms = table->atomCount();
            Selection sub(s.size());
            for (auto i=table->begin(), e=table->end(); i!=e; ++i) {
                if (std::binary_search(ids.begin(), ids.end(), i->param())) {
                    for (Id j=0; j<natoms; j++) {
                        sub[i->atom(j)] = 1;
                    }
                }
            }
            s.intersect(sub);
        }
        void dump(std::ostream& out) const {}
    };
}

StringPredicate* VMD::find_strfctn(const char* s) {
    StringPredicate* p;
    if      (!strcmp(s,"smarts")) p = new SmartsPredicate(sys);
    else if (!strcmp(s,"paramtype")) p = new ParamTypePredicate(sys);
    else    return NULL;
    predicates.push_back(PredicatePtr(p));
    return p;
}

Predicate* VMD::make_strfctn(StringPredicate* p, StringList* targets) {
    std::vector<std::string> values(targets->values.begin(), targets->values.end());
    delete targets;
    for (unsigned i=0; i<values.size(); i++) p->add(values[i]);
    return p;
}

Predicate* VMD::find_macro(const char* s) {
    unsigned i,n = sizeof(builtin_macros)/sizeof(builtin_macros[0]);
    const char* macro = NULL;
    for (i=0; i<n; i++) {
        if (!strcmp(s, builtin_macros[i].name)) {
            macro = builtin_macros[i].text;
            break;
        }
    }
    if (!macro) return NULL;
    return add(VMD(sys).parse(macro));
}

bool VMD::find_function(const char* s) {
    /* FIXME: must agree with atomsel::unary_expression */
    return !strcmp(s,"sqr") ||
           !strcmp(s,"sqrt") ||
           !strcmp(s,"abs");
}

Predicate* VMD::make_not(Predicate* sub) {
    return add(not_predicate(sub->shared_from_this()));
}

Predicate* VMD::make_and(Predicate* lhs, Predicate* rhs) {
    return add(and_predicate(lhs->shared_from_this(), rhs->shared_from_this()));
}

Predicate* VMD::make_or(Predicate* lhs, Predicate* rhs) {
    return add(or_predicate(lhs->shared_from_this(), rhs->shared_from_this()));
}

Predicate* VMD::make_single(Keyword* key) {
    return add(boolean_predicate(key->shared_from_this()));
}

Predicate* VMD::make_key(Keyword* key, TargetList* targets) {
    return add(keyword_predicate(key->shared_from_this(), targets));
}

Predicate* VMD::make_within(double r, Predicate* sub) {
    return add(within_predicate(sys, r, sub->shared_from_this()));
}
Predicate* VMD::make_pbwithin(double r, Predicate* sub) {
    return add(pbwithin_predicate(sys, r, sub->shared_from_this()));
}
Predicate* VMD::make_exwithin(double r, Predicate* sub) {
    return add(exwithin_predicate(sys, r, sub->shared_from_this()));
}
Predicate* VMD::make_withinbonds(int n, Predicate* sub) {
    return add(withinbonds_predicate(sys, n, sub->shared_from_this()));
}

Predicate* VMD::make_same(Keyword* k, Predicate* sub) {
    return add(same_keyword_predicate(k->shared_from_this(),
                                      full_selection(sys),
                                      sub->shared_from_this()));
}

Predicate* VMD::make_nearest(int r, Predicate* sub) {
    return add(k_nearest_predicate(sys, r, sub->shared_from_this()));
}

Predicate* VMD::make_pbnearest(int r, Predicate* sub) {
    return add(k_pbnearest_predicate(sys, r, sub->shared_from_this()));
}

Expression* VMD::make_exp(int v) {
    return add(literal_expression(v));
}

Expression* VMD::make_exp(double v) {
    return add(literal_expression(v));
}

Expression* VMD::make_exp(Keyword* k) {
    return add(keyword_expression(k->shared_from_this()));
}

Expression* VMD::make_unaexp(const char* op, Expression* s) {
    return add(unary_expression(op, s->shared_from_this()));
}

Expression* VMD::make_binexp(const char* op, Expression* lhs, Expression* rhs) {
    return add(binary_expression(op, lhs->shared_from_this(),
                                     rhs->shared_from_this()));
}

Predicate* VMD::make_compare(int cmp, Expression* lhs, Expression* rhs) {
    return add(relation_predicate(static_cast<RelationType>(cmp), 
                                  lhs->shared_from_this(),
                                  rhs->shared_from_this()));
}
