#ifndef desres_msys_elements_hxx
#define desres_msys_elements_hxx

namespace desres { namespace msys {

    int GuessAtomicNumber( double mass );
    const char* AbbreviationForElement(int anum);
    int ElementForAbbreviation(const char* abbr);

    double RadiusForElement(int anum);

    struct ChemData {
        float eneg;    // Electronegativities (Allen Scale)
        int nValence;  // Number of valence electrons present in atom
        int maxOct;    // max electrons for 'satisfied' octet (hypervalent>8)
        int maxFree;   // max free (non-bonding) electron pairs allowed...

        ChemData() 
        : eneg(), nValence(), maxOct(), maxFree() 
        {}

        ChemData(float const& e, int const& v, int const& o, int const& f)
        : eneg(e), nValence(v), maxOct(o), maxFree(f) 
        {}

        bool nodata() const { return eneg==0; }
    };

    /* return ChemData for the given atomic number, or invalid ChemData
     * (nodata()==true) if not available */
    ChemData const& DataForElement(int anum);
}}

#endif
