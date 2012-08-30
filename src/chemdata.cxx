#include "elements.hxx"

using namespace desres::msys;

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


