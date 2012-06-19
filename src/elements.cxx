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
