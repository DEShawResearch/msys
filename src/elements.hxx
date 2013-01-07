#ifndef desres_msys_elements_hxx
#define desres_msys_elements_hxx

namespace desres { namespace msys {

    int GuessAtomicNumber( double mass );
    const char* AbbreviationForElement(int anum);
    double MassForElement(int anum);
    int ElementForAbbreviation(const char* abbr);

    double RadiusForElement(int anum);

    struct ChemData {
        float eneg;    // Electronegativities (Allen Scale)
        unsigned nValence;  // Number of valence electrons present in atom
        unsigned maxOct;    // max electrons for 'satisfied' octet (hypervalent>8)
        unsigned maxFree;   // max free (non-bonding) electron pairs allowed...
        unsigned maxCoord;  // max coordinated atoms (prevent over-bonding of 1-2 row elements)

        ChemData() 
        : eneg(), nValence(), maxOct(), maxFree(), maxCoord()
        {}

        ChemData(float const& e, unsigned const& v, unsigned const& o, 
                 unsigned const& f, unsigned const& c)
        : eneg(e), nValence(v), maxOct(o), maxFree(f), maxCoord(c)
        {}

        bool nodata() const { return eneg==0; }
    };

    /* return ChemData for the given atomic number, or invalid ChemData
     * (nodata()==true) if not available */
    ChemData const& DataForElement(int anum);
}}

#endif
