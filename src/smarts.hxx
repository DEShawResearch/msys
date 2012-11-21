#ifndef desres_msys_smarts_hxx
#define desres_msys_smarts_hxx

#include "system.hxx"

namespace desres { namespace msys {

    /* The smarts pattern implementation has lots of hairy boost::spirit
     * stuff that we'd prefer not to expose.
     */
    class SmartsPatternImpl;
    typedef boost::shared_ptr<SmartsPatternImpl> SmartsPatternImplPtr;

    /* The SmartsPattern class holds a parsed SMARTS pattern and can
     * efficiently find matches of that pattern in a given query 
     * structure.
     *
     * This class is not a shared pointer; it can be efficiently copied.
     */
    class SmartsPattern {
        std::string             _pattern;
        SmartsPatternImplPtr    _impl;
        std::string             _warnings;

    public:
        explicit SmartsPattern(std::string const& pattern);
        Id atomCount() const;
        std::string const& pattern() const { return _pattern; }

        /* any warning messages encountered during processing of the pattern */
        std::string const& warnings() const { return _warnings; }

        /* add atom and properties to the given system; these must be
         * present before calling match().  
         *
         * Assigns bond orders, formal charges, aromaticity, ring sizes, ring 
         * counts, hydrogen counts, bond degrees.
         */
        static void Annotate(SystemPtr sys, IdList const& atoms);

        /* Find matches of the SMARTS pattern that start with the given
         * set of atoms; usually sys->atoms().  System must have atom
         * and bond properties populated by GenerateAtomAndBondProperties.
         * Will give duplicate matches for '*~*' in forward and reverse
         * ordering, etc. */
         MultiIdList findMatches(SystemPtr sys, IdList const& starts) const;
    };

}}

#endif
