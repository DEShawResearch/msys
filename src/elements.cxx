#include "elements.hxx"
#include <algorithm>
#include <string.h>

static const char *ename[] = {
    "X",  "H",  "He", "Li", "Be", "B",  "C",  "N",  "O",  "F",  "Ne",
    "Na", "Mg", "Al", "Si", "P" , "S",  "Cl", "Ar", "K",  "Ca", "Sc",
    "Ti", "V",  "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge",
    "As", "Se", "Br", "Kr", "Rb", "Sr", "Y",  "Zr", "Nb", "Mo", "Tc",
    "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I",  "Xe",
    "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb",
    "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W",  "Re", "Os",
    "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr",
    "Ra", "Ac", "Th", "Pa", "U",  "Np", "Pu", "Am", "Cm", "Bk", "Cf",
    "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt",
    "Ds", "Rg"
};

/* atomic weights do not increase monotonically!  I'm writing down
 * the masses in order of atomic number, and I'll then assign
 * atomic number and then sort by mass. */
namespace {
    struct Element {
        double  mass;
        int     anum;
        Element() {}
        explicit Element(double m) : mass(m), anum(0) {}
        Element(double m, int a) : mass(m), anum(a) {}
        bool operator<(Element const& e) const { return mass<e.mass; }
    };

    static Element elems[] = {
        {0.00000,    0},
        {1.00794,    0},
        {4.00260,    0},
        {6.941,      0},
        {9.012182,   0},
        {10.811,     0},

        {12.0107,    0},
        {14.0067,    0},
        {15.9994,    0},
        {18.9984032, 0},
        {20.1797,    0},

        {22.989770,  0},
        {24.3050,    0},
        {26.981538,  0},
        {28.0855,    0},
        {30.973761,  0},

        {32.065,     0},
        {35.453,     0},
        {39.948,     0},
        {39.0983,    0},
        {40.078,     0},
        {44.955910,  0},

        {47.867,     0},
        {50.9415,    0},
        {51.9961,    0},
        {54.938049,  0},
        {55.845,     0},
        {58.9332,    0},

        {58.6934,    0},
        {63.546,     0},
        {65.409,     0},
        {69.723,     0},
        {72.64,      0},
        {74.92160,   0},

        {78.96,  0},
        {79.904, 0},
        {83.798, 0},
        {85.4678,    0},
        {87.62,  0},
        {88.90585,   0},

        {91.224, 0},
        {92.90638,   0},
        {95.94,  0},
        {98.0,   0},
        {101.07, 0},
        {102.90550,  0},

        {106.42, 0},
        {107.8682,   0},
        {112.411,    0},
        {114.818,    0},
        {118.710,    0},
        {121.760,    0},

        {127.60, 0},
        {126.90447,  0},
        {131.293,    0},
        {132.90545,  0},
        {137.327,    0},

        {138.9055,   0},
        {140.116,    0},
        {140.90765,  0},
        {144.24, 0},
        {145.0,  0},
        {150.36, 0},

        {151.964,    0},
        {157.25, 0},
        {158.92534,  0},
        {162.500,    0},
        {164.93032,  0},

        {167.259,    0},
        {168.93421,  0},
        {173.04, 0},
        {174.967,    0},
        {178.49, 0},
        {180.9479,   0},

        {183.84, 0},
        {186.207,    0},
        {190.23, 0},
        {192.217,    0},
        {195.078,    0},
        {196.96655,  0},

        {200.59, 0},
        {204.3833,   0},
        {207.2,  0},
        {208.98038,  0},
        {209.0,  0},
        {210.0,  0},
        {222.0,  0},

        {223.0,  0},
        {226.0,  0},
        {227.0,  0},
        {232.0381,   0},
        {231.03588,  0},
        {238.02891,  0},

        {237.0,  0},
        {244.0,  0},
        {243.0,  0},
        {247.0,  0},
        {247.0,  0},
        {251.0,  0},
        {252.0,  0},
        {257.0,  0},

        {258.0,  0},
        {259.0,  0},
        {262.0,  0},
        {261.0,  0},
        {262.0,  0},
        {266.0,  0},
        {264.0,  0},
        {269.0,  0},

        {268.0,  0},
        {271.0,  0},
        {272.0,  0}
    };

    static const unsigned nelems = sizeof(elems)/sizeof(elems[0]);

    struct eleminit {
        eleminit() {
            for (unsigned i=0; i<nelems; i++) elems[i].anum=i;
            std::sort(elems, elems+nelems);
        }
    } ini;
}

int desres::msys::GuessAtomicNumber( double mass ) {
    const Element* begin=elems, *end = elems+nelems;
    /* first element which is not less than mass */
    const Element* rhs = std::lower_bound(begin, end, Element(mass));
    /* greater than last value */
    if (rhs==end) return 0;
    /* zero or negative mass? */
    if (rhs==begin) return 0;
    /* take the closest */
    const Element* lhs = rhs-1;
    double ldiff = mass-lhs->mass;
    double rdiff = rhs->mass-mass;
    return ldiff < rdiff ? lhs->anum : rhs->anum;
}

const char* desres::msys::AbbreviationForElement(int anum) {
    int n = sizeof(ename)/sizeof(ename[0]);
    if (anum<0 || anum>=n) return "";
    return ename[anum];
}
