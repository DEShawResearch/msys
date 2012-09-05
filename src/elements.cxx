#include "elements.hxx"
#include <algorithm>
#include <string.h>
#include <string>
#include <boost/algorithm/string.hpp>

using namespace desres::msys;

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

static double eradius[] = {
    0.00, 1.10, 1.40, 1.81, 1.53, 1.92, 1.70, 1.55, 1.52, 1.47, 1.54,
    1.81, 1.73, 1.84, 2.10, 1.80, 1.80, 1.75, 1.88, 2.75, 2.31, 2,
    2,    2,    2,    2,    2,    2,    2,    2,    2,    1.87, 2.11,
    1.85, 1.90, 1.83, 2.02, 3.03, 2.49, 2,    2,    2,    2,    2,
    2,    2,    2,    2,    2,    1.93, 2.17, 2.06, 2.06, 1.98, 2.16,
    3.43, 2.68, 2,    2,    2,    2,    2,    2,    2,    2,    2,
    2,    2,    2,    2,    2,    2,    2,    2,    2,    2,    2,
    2,    2,    2,    2,    1.96, 2.02, 2.07, 1.97, 2.02, 2.20, 3.48,
    2.83, 2,    2,    2,    2,    2,    2,    2,    2,    2,    2,
    2,    2
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
        Element(0.00000,    0),
        Element(1.00794,    0),
        Element(4.00260,    0),
        Element(6.941,      0),
        Element(9.012182,   0),
        Element(10.811,     0),

        Element(12.0107,    0),
        Element(14.0067,    0),
        Element(15.9994,    0),
        Element(18.9984032, 0),
        Element(20.1797,    0),

        Element(22.989770,  0),
        Element(24.3050,    0),
        Element(26.981538,  0),
        Element(28.0855,    0),
        Element(30.973761,  0),

        Element(32.065,     0),
        Element(35.453,     0),
        Element(39.948,     0),
        Element(39.0983,    0),
        Element(40.078,     0),
        Element(44.955910,  0),

        Element(47.867,     0),
        Element(50.9415,    0),
        Element(51.9961,    0),
        Element(54.938049,  0),
        Element(55.845,     0),
        Element(58.9332,    0),

        Element(58.6934,    0),
        Element(63.546,     0),
        Element(65.409,     0),
        Element(69.723,     0),
        Element(72.64,      0),
        Element(74.92160,   0),

        Element(78.96,  0),
        Element(79.904, 0),
        Element(83.798, 0),
        Element(85.4678,    0),
        Element(87.62,  0),
        Element(88.90585,   0),

        Element(91.224, 0),
        Element(92.90638,   0),
        Element(95.94,  0),
        Element(98.0,   0),
        Element(101.07, 0),
        Element(102.90550,  0),

        Element(106.42, 0),
        Element(107.8682,   0),
        Element(112.411,    0),
        Element(114.818,    0),
        Element(118.710,    0),
        Element(121.760,    0),

        Element(127.60, 0),
        Element(126.90447,  0),
        Element(131.293,    0),
        Element(132.90545,  0),
        Element(137.327,    0),

        Element(138.9055,   0),
        Element(140.116,    0),
        Element(140.90765,  0),
        Element(144.24, 0),
        Element(145.0,  0),
        Element(150.36, 0),

        Element(151.964,    0),
        Element(157.25, 0),
        Element(158.92534,  0),
        Element(162.500,    0),
        Element(164.93032,  0),

        Element(167.259,    0),
        Element(168.93421,  0),
        Element(173.04, 0),
        Element(174.967,    0),
        Element(178.49, 0),
        Element(180.9479,   0),

        Element(183.84, 0),
        Element(186.207,    0),
        Element(190.23, 0),
        Element(192.217,    0),
        Element(195.078,    0),
        Element(196.96655,  0),

        Element(200.59, 0),
        Element(204.3833,   0),
        Element(207.2,  0),
        Element(208.98038,  0),
        Element(209.0,  0),
        Element(210.0,  0),
        Element(222.0,  0),

        Element(223.0,  0),
        Element(226.0,  0),
        Element(227.0,  0),
        Element(232.0381,   0),
        Element(231.03588,  0),
        Element(238.02891,  0),

        Element(237.0,  0),
        Element(244.0,  0),
        Element(243.0,  0),
        Element(247.0,  0),
        Element(247.0,  0),
        Element(251.0,  0),
        Element(252.0,  0),
        Element(257.0,  0),

        Element(258.0,  0),
        Element(259.0,  0),
        Element(262.0,  0),
        Element(261.0,  0),
        Element(262.0,  0),
        Element(266.0,  0),
        Element(264.0,  0),
        Element(269.0,  0),

        Element(268.0,  0),
        Element(271.0,  0),
        Element(272.0,  0)
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

double desres::msys::RadiusForElement(int anum) {
    int n = sizeof(eradius)/sizeof(eradius[0]);
    if (anum<0) return 0;
    if (anum>=n) return 2.0;
    return eradius[anum];
}

int desres::msys::ElementForAbbreviation(const char* abbr) {
    std::string src(abbr);
    boost::to_upper(src);
    for (unsigned i=1; i<nelems; i++) {
        std::string ref(ename[i]);
        boost::to_upper(ref);
        if (ref==src) return i;
    }
    return 0;
}

namespace {
    ChemData NODATA;

    const ChemData atomInfo[] = { 
        NODATA,                                 // PlaceHolder
     ChemData(2.300,   1,       2,       0 ),        // H   1
     ChemData(4.160,   2,       2,       0 ),        // He  2
        // Period 2
     ChemData(0.912,   1,       8,       0 ),        // Li  3
     ChemData(1.576,   2,       8,       0 ),        // Be  4
     ChemData(2.051,   3,       8,       0 ),        // B   5
     ChemData(2.544,   4,       8,       1 ),        // C   6
     ChemData(3.066,   5,       8,       2 ),        // N   7
     ChemData(3.610,   6,       8,       3 ),        // O   8
     ChemData(4.193,   7,       8,       4 ),        // F   9
     ChemData(4.789,   8,       8,       4 ),        // Ne  10
        // Period 3
     ChemData(0.869,   1,       8,       0 ),        // Na  11
     ChemData(1.293,   2,       8,       0 ),        // Mg  12
     ChemData(1.613,   3,       8,       0 ),        // Al  13
     ChemData(1.916,   4,       8,       1 ),        // Si  14
     ChemData(2.253,   5,       10,      2 ),        // P   15
     ChemData(2.589,   6,       12,      3 ),        // S   16
     ChemData(2.869,   7,       8,       4 ),        // Cl  17
     ChemData(3.242,   8,       8,       4 ),        // Ar  18
        // Period 4 - Groups 1 & 2
     ChemData( 0.734,  1,       8,       0 ),        // K   19
     ChemData( 1.034,  2,       8,       0 ),        // Ca  20
        // Period 4 - Transition Metals (most stable or lowest charge state is used)
     ChemData(1.190,   3,       8,       0 ),        // Sc  21 (2 rare, 3 common)
     ChemData(1.380,   4,       8,       0 ),        // Ti  22 (2-3 rare, 4 common)
     ChemData(1.530,   2,       8,       0 ),        // V   23 (1 rare, 2-5 common)
     ChemData(1.650,   3,       8,       0 ),        // Cr  24 (1,4-5 rare, 2-3,6 common, 3 most stable) 
     ChemData(1.750,   2,       8,       0 ),        // Mn  25 (1,5 rare, 2-4,6-7 common, 2 most stable)
     ChemData(1.800,   2,       8,       0 ),        // Fe  26 (1,4-6 rare, 2-3 common)
     ChemData(1.840,   2,       8,       0 ),        // Co  27 (1,4-6 rare, 2-3 common)
     ChemData(1.880,   2,       8,       0 ),        // Ni  28 (1,3-4 rare, 2 common)
     ChemData(1.850,   2,       8,       0 ),        // Cu  29 (1 rare, 2 common)
     ChemData(1.590,   2,       8,       0 ),        // Zn  30 (2 common)
        // Period 4 - Groups 13-18
        NODATA, NODATA,                         // Ga - Ge    31 - 32
     ChemData( 2.211,  5,       10,      2 ),        // As  33
     ChemData( 2.434,  6,       12,      3 ),        // Se  34
     ChemData( 2.685,  7,       8,       4 ),        // Br  35
     ChemData( 2.966,  8,       8,       4 ),        // Kr  36
        // Period 5 - Groups 1 & 2
     ChemData( 0.706,  1,       8,       0 ),        // Rb  37
     ChemData( 0.963,  2,       8,       0 ),        // Sr  38
        // Period 5 - Transition Metals
        // FIXME: should handle biocatalytic metals (Mo) 
        NODATA, NODATA, NODATA, NODATA, NODATA, //  Y - Tc  39-43
        NODATA, NODATA, NODATA, NODATA, NODATA, // Ru - Cd  44-48
        // Period 5 - Groups 13-18
        NODATA, NODATA, NODATA, NODATA,         // In - Te  49-52
     ChemData( 2.359,  7,       8,       4 ),        // I   53
     ChemData( 2.582,  8,       8,       4 ),        // Kr  54          
        // Period 6 - Groups 1 & 2
     ChemData( 0.659,  1,       8,       0 ),         // Cs  55
     ChemData( 0.881,  2,       8,       0 ),         // Ba  56
        // Period 6 - Lanthanides (all form M+3, some form M+2 or M+4)
        NODATA, NODATA, NODATA, NODATA, NODATA, NODATA, NODATA, // La - Eu  57 - 63
        NODATA, NODATA, NODATA, NODATA, NODATA, NODATA, NODATA, // Gd - Yb  64 - 70
        // Period 6 - Transition Metals
        NODATA, NODATA, NODATA, NODATA, NODATA, // Lu - Re 71-75
        NODATA, NODATA, NODATA, NODATA, NODATA, // Os - Hg 76-80
        // Period 6 - Groups 13-18
        NODATA, NODATA, NODATA, NODATA, NODATA, // Ti - At 81-85
     ChemData( 2.60,   8,       8,       4 ),        // Rn  86     
        // Period 7 - Groups 1 & 2
     ChemData( 0.659,  1,       8,       0 ),        // Fr  87
     ChemData( 0.881,  2,       8,       0 )         // Ra  88
    };

    const int maxData = sizeof(atomInfo)/sizeof(atomInfo[0]);
}

ChemData const& desres::msys::DataForElement(int anum) {
    if (anum<0 || anum>=maxData) return NODATA;
    return atomInfo[anum];
}


