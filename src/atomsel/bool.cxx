#include "token.hxx"
#include <unordered_map>

using namespace desres::msys;
using namespace desres::msys::atomsel;

typedef void (*eval_t)(Selection& s, System* mol);

static void eval_all(Selection& s, System* mol) {
    s.fill();
}

static void eval_none(Selection& s, System* mol) {
    s.clear();
}

static void eval_hydrogen(Selection& s, System* mol) {
    for (Id i=0, n=s.size(); i<n; i++) {
        s[i]=s[i] && mol->atomFAST(i).atomic_number==1;
    }
}
static void eval_oxygen(Selection& s, System* mol) {
    for (Id i=0, n=s.size(); i<n; i++) {
        s[i]=s[i] && mol->atomFAST(i).atomic_number==8;
    }
}

static void eval_noh(Selection& s, System* mol) {
    for (Id i=0, n=s.size(); i<n; i++) {
        s[i]=s[i] && mol->atomFAST(i).atomic_number!=1;
    }
}

static void eval_backbone(Selection& s, System* mol) {
    for (Id i=0, n=s.size(); i<n; i++) {
        s[i]=s[i] && (mol->atomFAST(i).type == AtomProBack ||
                      mol->atomFAST(i).type == AtomNucBack);
    }
}

static void eval_sidechain(Selection& s, System* mol) {
    for (Id i=0, n=s.size(); i<n; i++) {
        s[i]=s[i] && mol->atomFAST(i).type == AtomProSide;
    }
}

static void eval_protein(Selection& s, System* mol) {
    for (Id i=0, n=s.size(); i<n; i++) {
        s[i]=s[i] && mol->residueFAST(mol->atomFAST(i).residue).type==ResidueProtein;
    }
}

static void eval_nucleic(Selection& s, System* mol) {
    for (Id i=0, n=s.size(); i<n; i++) {
        s[i]=s[i] && mol->residueFAST(mol->atomFAST(i).residue).type==ResidueNucleic;
    }
}

static void eval_water(Selection& s, System* mol) {
    for (Id i=0, n=s.size(); i<n; i++) {
        s[i]=s[i] && mol->residueFAST(mol->atomFAST(i).residue).type==ResidueWater;
    }
}

static const std::unordered_map<std::string,eval_t> map = {
    {"all", eval_all},
    {"none", eval_none},
    {"hydrogen", eval_hydrogen},
    {"oxygen", eval_oxygen},
    {"noh", eval_noh},
    {"backbone", eval_backbone},
    {"sidechain", eval_sidechain},
    {"protein", eval_protein},
    {"nucleic", eval_nucleic},
    {"water", eval_water},
};

void BoolPredicate::eval(Selection& s) {
    auto p = map.find(name);
    if (p==map.end()) {
        MSYS_FAIL("Unrecognized boolean '" << name << "'");
    }
    p->second(s, mol);
}
