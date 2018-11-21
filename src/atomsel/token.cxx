#include "token.hxx"
#include <unordered_map>
#include <cmath>

using namespace desres::msys::atomsel;

static std::unordered_map<std::string,int> reserved = {
    {"of", OF},
    {"same", SAME},
    {"as", AS},
    {"to", TO},
    {"and", AND},
    {"or", OR},
    {"not", NOT},
    {"within", WITHIN},
    {"exwithin", EXWITHIN},
    {"pbwithin", PBWITHIN},
    {"withinbonds", WITHINBONDS},
    {"nearest", NEAREST},
    {"pbnearest", PBNEAREST},
};

static double sqr(double x) { return x*x; }

typedef double (*func_t)(double);
static const std::unordered_map<std::string,func_t> funcs = {
    {"sqr", sqr},
    {"sqrt", std::sqrt},
    {"abs", std::fabs},
};

static std::unordered_map<std::string,std::string> macros = {
    {"at","resname ADE A THY T"},
    {"acidic","resname ASP GLU"},
    {"cyclic","resname HIS PHE PRO TRP TYR"},
    {"acyclic","protein and not cyclic"},
    {"aliphatic","resname ALA GLY ILE LEU VAL"},
    {"alpha","protein and name CA"},
    {"amino","protein"},
    {"aromatic","resname HIS PHE TRP TYR"},
    {"basic","resname ARG HIS LYS HSP"},
    {"bonded","degree > 0"},
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
    {"legacy_ion","resname AL BA CA Ca CAL CD CES CLA CL Cl 'Cl-' CO CS CU Cu CU1 CUA HG IN IOD K 'K+' MG MN3 MO3 MO4 MO5 MO6 NA Na NAW OC7 PB POT PT RB SOD TB TL WO4 YB ZN ZN1 ZN2"},
    // there are too many possible atomic numbers to list.  both amber and charmm have parameters
    // for dozens of monoatomic species and even a few diatomics.  So just exclude organic
    // atoms which are probably never intended to be bare ions, and the noble gases.
    {"ion", "degree 0 and not atomicnumber 0 1 2 5 6 7 8 10 18 36 54 86" },
    {"ions","ion"},
    {"sugar","resname AGLC"},
    {"solvent","not (protein or sugar or nucleic or lipid)"},
    {"carbon","atomicnumber 6"},
    {"nitrogen","atomicnumber 7"},
    {"sulfur","atomicnumber 16"},
    {"heme","resname HEM HEME"},
};

static const char syntax[] = " '\"()<=>+-%*/^\t\n";

int Tokenizer::next(desres::msys::atomsel::Token* t, desres::msys::System* mol) {
    // check for eof
    if (!s[loc]) return 0;

    // skip space
    while (s[loc] && isspace(s[loc])) ++loc;
    if (!s[loc]) return 0;  // allow trailing space

    // skip ahead until we find syntax characters
    int wordlen = strcspn(s+loc, syntax);

    // did we skip past a possible token?
    if (wordlen>0) {
        std::string text(s+loc,s+loc+wordlen);
        t->sval = s+loc;
        t->slen = wordlen;
        loc += wordlen;

        // allow single quote in the word (5', C4', 3'A)
        if (s[loc]=='\'') {
            // accept additional characters until we hit either end of
            // string or syntax.  If the next syntax character is yet
            // another single quote, it's an error.
            auto extra = strcspn(s+loc+1, syntax)+1;
            t->slen += extra;
            if (s[t->slen]=='\'') {
                MSYS_FAIL("Malformed string value \"" << t->str() << "\"");
            }
            loc += extra;
            return VAL;
        }

        // check for reserved words
        auto p = reserved.find(text);
        if (p!=reserved.end()) return p->second;

        // check for functions.
        {
            auto p = funcs.find(text);
            if (p!=funcs.end() && s[loc]=='(') {
                t->func = p->second;
                return FUNC;
            }
        }

        // check for macros
        {
            auto p = macros.find(text);
            if (p!=macros.end()) {
                t->sval = p->second.data();
                t->slen = p->second.size();
                return MACRO;
            }
        }

        // check for int
        char* end;
        errno = 0;
        t->ival = strtol(text.data(), &end, 10);
        if (errno==ERANGE) {
            //MSYS_FAIL("Integer value " << s << " would overflow");
            return -1;
        }
        if (*end==0) {
            return INT;
        }

        // check for float
        errno = 0;
        t->fval = strtod(text.data(), &end);
        if (errno==ERANGE) {
            //MSYS_FAIL("float value " << s << " would overflow");
            return -1;
        }
        if (*end==0) {
            return FLT;
        }

        // check for keyword
        if (is_keyword(text, mol)) {
            return KEY;
        }

        // regular string
        return VAL;
    }

    // check for quoted words.  Double quote signifies regex; single quote
    // for a literal string.
    if (s[loc]=='"' || s[loc]=='\'') {
        char quote = s[loc];
        ++loc;
        t->sval = s+loc;
        // search for matching unescaped close quote
        while (s[loc]) {
            if (s[loc]==quote && s[loc-1] != '\\') break;
            ++loc;
        }
        if (s[loc]==quote) {
            t->slen = (s+loc)-t->sval;
            ++loc;
            return quote=='"' ? REGEX : VAL;
        }
        // unterminated quote
        return -1;
    }

    // just have syntax to return
    if (s[loc]=='(') { ++loc; return LPAREN; }
    if (s[loc]==')') { ++loc; return RPAREN; }

    // comparison operators
    if (s[loc]=='<') {
        if (s[loc+1]=='=') {
            t->ival = Token::LE;
            ++loc;
        } else {
            t->ival = Token::LT;
        }
        ++loc;
        return CMP;
    }
    if (s[loc]=='>') {
        if (s[loc+1]=='=') {
            t->ival = Token::GE;
            ++loc;
        } else {
            t->ival = Token::GT;
        }
        ++loc;
        return CMP;
    }
    if (s[loc]=='=' && s[loc+1]=='=') {
        loc += 2;
        t->ival = Token::EQ;
        return CMP;
    }
    if (s[loc]=='!' && s[loc+1]=='=') {
        loc += 2;
        t->ival = Token::NE;
        return CMP;
    }

    // plus and minus have special handling.  If the character is at the
    // start of the selection, or if syntax precedes the character and
    // a non-syntax character or paren follows, they're unary; otherwise binary.
    if (s[loc]=='-') {
        if (loc==0 || (strchr(syntax,s[loc-1]) && 
                      (s[loc+1]=='(' || !strchr(syntax,s[loc+1])))) {
            ++loc;
            return NEG;
        }
        ++loc;
        return SUB;
    }
    if (s[loc]=='+') {
        if (loc==0 || (strchr(syntax,s[loc-1]) && (s[loc+1]=='.' || isdigit(s[loc+1])))) {
            ++loc;
            return PLUS;
        }
        ++loc;
        return ADD;
    }
    // other binary operators
    if (s[loc]=='*') {
        ++loc;
        if (s[loc+1]=='*') {
            ++loc;
            return EXP;
        }
        return MUL;
    }
    if (s[loc]=='/') {
        ++loc;
        return DIV;
    }
    if (s[loc]=='^') {
        ++loc;
        return EXP;
    }
    if (s[loc]=='%') {
        ++loc;
        return MOD;
    }

    // unrecognized syntax
    return -1;
}

