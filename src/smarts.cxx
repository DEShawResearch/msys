#include "smarts.hxx"
#include <boost/fusion/adapted/struct/adapt_struct.hpp>
#include <boost/fusion/include/adapt_struct.hpp>
#include <boost/spirit/include/phoenix.hpp>
#include <boost/spirit/include/qi.hpp>
#include <boost/tuple/tuple.hpp>
#include <stack>

namespace qi = boost::spirit::qi;
namespace ascii = boost::spirit::ascii;
using namespace desres::msys;

/* Return types for boost::spirit::qi parser */
namespace { 
    namespace bf = boost::fusion;
    struct smarts_pattern_;
    enum AROMATICITY { AROM, ALIPH, NONE };
    typedef boost::variant<std::pair<int, AROMATICITY>, char> raw_atom_;
    typedef std::pair<int, AROMATICITY> element_;
    typedef bf::vector2<char, unsigned> atomic_number_;
    typedef bf::vector2<char, boost::optional<unsigned> >
        optional_numeric_property_;
    typedef boost::variant<bf::vector2<char, unsigned>,
            std::vector<char> > charge_;
    typedef boost::variant<boost::recursive_wrapper<smarts_pattern_>,
            element_, atomic_number_, optional_numeric_property_,
            charge_, SmartsPatternImplPtr> atom_spec_;
    typedef bf::vector2<std::string, atom_spec_> not_atom_spec_;
    typedef std::vector<std::vector<not_atom_spec_> > and_atom_spec_;
    typedef std::vector<and_atom_spec_> or_atom_spec_;
    typedef std::vector<or_atom_spec_> atom_expression_;
    typedef bf::vector2<char, boost::optional<charge_> >
        hydrogen_expression_;
    typedef boost::variant<raw_atom_, hydrogen_expression_,
            atom_expression_> atom_;
    typedef std::string bond_spec_;
    typedef std::vector<std::vector<bond_spec_> > and_bond_spec_;
    typedef std::vector<and_bond_spec_> or_bond_spec_;
    typedef boost::optional<std::vector<or_bond_spec_> > bond_expression_;
    typedef bf::vector2<bond_expression_,
            boost::recursive_wrapper<smarts_pattern_> > bond_;
    typedef boost::variant<char, bf::vector2<char, char> > closure_ind_;
    typedef std::vector<bf::vector2<bond_expression_, int> > closure_;
    struct smarts_pattern_ {
        atom_ atom;
        closure_ closure;
        std::vector<bond_> bonds;
        boost::optional<bond_> last_bond;
    };
}

namespace desres { namespace msys {

    class SmartsPatternImpl {

        /* All atom expressions, in their original order of specification */
        std::vector<atom_> _atoms;

        /* (atom_index, bond_expression, atom_index) for all bond
         * expressions, in their original order of specification, with
         * closure bonds appearing at the location of their second
         * specification. Closure bonds have second atom index less than
         * first atom index; all other bonds have second atom index
         * corresponding to a new unexplored atom and greater than first
         * atom index. */
        std::vector<boost::tuple<unsigned,
            bond_expression_, unsigned> >
            _bonds;

        /* Helper function to construct _atoms and _bonds lists from
         * smarts_pattern_ structure returned by parser */
        bool convertSmarts(const smarts_pattern_&
                pattern, std::map<int,
                std::pair<bond_expression_,
                unsigned> >& closure_map);

        /* construct only through ::create(). */
        SmartsPatternImpl(std::string const& pat, std::ostream& log);

    public:
        static SmartsPatternImplPtr create(std::string const& pat,
                              std::ostream& log)  {
            return SmartsPatternImplPtr(new SmartsPatternImpl(pat, log));
        }

        Id atomCount() const { return _atoms.size(); }

        /* Helper function to match a SMARTS pattern starting at a given
         * atom and append matches to a given list. Has option of returning
         * after finding a single match. Returns true if any match is found,
         * false otherwise. */
        bool matchSmartsPattern(AnnotatedSystemPtr sys, Id atom,
                MultiIdList& matches, bool match_single=false) const;
    };

}}

/************************ Symbol tables for elements **************************/

struct raw_single_sym_ : qi::symbols<char, std::pair<int, AROMATICITY> > {
    /* Single character, can be specified outside [] */
    raw_single_sym_() {
        add
            ("*", std::make_pair(0, NONE))
            ("A", std::make_pair(0, ALIPH))
            ("a", std::make_pair(0, AROM))
            ("B", std::make_pair(5, ALIPH))
            ("b", std::make_pair(5, AROM))
            ("C", std::make_pair(6, ALIPH))
            ("c", std::make_pair(6, AROM))
            ("N", std::make_pair(7, ALIPH))
            ("n", std::make_pair(7, AROM))
            ("O", std::make_pair(8, ALIPH))
            ("o", std::make_pair(8, AROM))
            ("F", std::make_pair(9, NONE))
            ("P", std::make_pair(15, ALIPH))
            ("p", std::make_pair(15, AROM))
            ("S", std::make_pair(16, ALIPH))
            ("s", std::make_pair(16, AROM))
            ("I", std::make_pair(53, NONE));
    }
} raw_single_sym;
struct raw_double_sym_ : qi::symbols<char, std::pair<int, AROMATICITY> > {
    /* Double character, can be specified outside [] */
    raw_double_sym_() {
        add
            ("Cl", std::make_pair(17, NONE))
            ("Se", std::make_pair(34, ALIPH))
            ("se", std::make_pair(34, AROM))
            ("Br", std::make_pair(35, NONE));
    }
} raw_double_sym;
struct single_sym_ : qi::symbols<char, std::pair<int, AROMATICITY> > {
    /* Single character, must be specified inside [] */
    single_sym_() {
        add
            ("K", std::make_pair(19, NONE))
            ("V", std::make_pair(23, NONE))
            ("Y", std::make_pair(39, NONE))
            ("W", std::make_pair(74, NONE))
            ("U", std::make_pair(92, NONE));
    }
} single_sym;
struct double_sym_ : qi::symbols<char, std::pair<int, AROMATICITY> >
{
    /* Double character, must be specified inside [] */
    double_sym_() {
        add
            ("He", std::make_pair(2, NONE))
            ("Li", std::make_pair(3, NONE))
            ("Be", std::make_pair(4, NONE))
            ("Ne", std::make_pair(10, NONE))
            ("Na", std::make_pair(11, NONE))
            ("Mg", std::make_pair(12, NONE))
            ("Al", std::make_pair(13, NONE))
            ("Si", std::make_pair(14, NONE))
            ("Ar", std::make_pair(18, NONE))
            ("Ca", std::make_pair(20, NONE))
            ("Sc", std::make_pair(21, NONE))
            ("Ti", std::make_pair(22, NONE))
            ("Cr", std::make_pair(24, NONE))
            ("Mn", std::make_pair(25, NONE))
            ("Fe", std::make_pair(26, NONE))
            ("Co", std::make_pair(27, NONE))
            ("Ni", std::make_pair(28, NONE))
            ("Cu", std::make_pair(29, NONE))
            ("Zn", std::make_pair(30, NONE))
            ("Ga", std::make_pair(31, NONE))
            ("Ge", std::make_pair(32, NONE))
            ("As", std::make_pair(33, ALIPH))
            ("as", std::make_pair(33, AROM))
            ("Kr", std::make_pair(36, NONE))
            ("Rb", std::make_pair(37, NONE))
            ("Sr", std::make_pair(38, NONE))
            ("Zr", std::make_pair(40, NONE))
            ("Nb", std::make_pair(41, NONE))
            ("Mo", std::make_pair(42, NONE))
            ("Tu", std::make_pair(43, NONE))
            ("Ru", std::make_pair(44, NONE))
            ("Rh", std::make_pair(45, NONE))
            ("Pd", std::make_pair(46, NONE))
            ("Ag", std::make_pair(47, NONE))
            ("Cd", std::make_pair(48, NONE))
            ("In", std::make_pair(49, NONE))
            ("Sn", std::make_pair(50, NONE))
            ("Sb", std::make_pair(51, NONE))
            ("Te", std::make_pair(52, NONE))
            ("Xe", std::make_pair(54, NONE))
            ("Cs", std::make_pair(55, NONE))
            ("Ba", std::make_pair(56, NONE))
            ("La", std::make_pair(57, NONE))
            ("Ce", std::make_pair(58, NONE))
            ("Pr", std::make_pair(59, NONE))
            ("Nd", std::make_pair(60, NONE))
            ("Pm", std::make_pair(61, NONE))
            ("Sm", std::make_pair(62, NONE))
            ("Eu", std::make_pair(63, NONE))
            ("Gd", std::make_pair(64, NONE))
            ("Tb", std::make_pair(65, NONE))
            ("Dy", std::make_pair(66, NONE))
            ("Ho", std::make_pair(67, NONE))
            ("Er", std::make_pair(68, NONE))
            ("Tm", std::make_pair(69, NONE))
            ("Yb", std::make_pair(70, NONE))
            ("Lu", std::make_pair(71, NONE))
            ("Hf", std::make_pair(72, NONE))
            ("Ta", std::make_pair(73, NONE))
            ("Re", std::make_pair(75, NONE))
            ("Os", std::make_pair(76, NONE))
            ("Ir", std::make_pair(77, NONE))
            ("Pt", std::make_pair(78, NONE))
            ("Au", std::make_pair(79, NONE))
            ("Hg", std::make_pair(80, NONE))
            ("Tl", std::make_pair(81, NONE))
            ("Pb", std::make_pair(82, NONE))
            ("Bi", std::make_pair(83, NONE))
            ("Po", std::make_pair(84, NONE))
            ("At", std::make_pair(85, NONE))
            ("Rn", std::make_pair(86, NONE))
            ("Fr", std::make_pair(87, NONE))
            ("Ra", std::make_pair(88, NONE))
            ("Ac", std::make_pair(89, NONE))
            ("Th", std::make_pair(90, NONE))
            ("Pa", std::make_pair(91, NONE))
            ("Np", std::make_pair(93, NONE))
            ("Pu", std::make_pair(94, NONE))
            ("Am", std::make_pair(95, NONE))
            ("Cm", std::make_pair(96, NONE))
            ("Bk", std::make_pair(97, NONE))
            ("Cf", std::make_pair(98, NONE))
            ("Es", std::make_pair(99, NONE))
            ("Fm", std::make_pair(100, NONE))
            ("Md", std::make_pair(101, NONE))
            ("No", std::make_pair(102, NONE))
            ("Lr", std::make_pair(103, NONE))
            ("Rf", std::make_pair(104, NONE))
            ("Db", std::make_pair(105, NONE))
            ("Sg", std::make_pair(106, NONE))
            ("Bh", std::make_pair(107, NONE))
            ("Hs", std::make_pair(108, NONE))
            ("Mt", std::make_pair(109, NONE))
            ("Ds", std::make_pair(110, NONE))
            ("Rg", std::make_pair(111, NONE));
    }
} double_sym;

/******************** boost::spirit::qi grammar and parser ********************/
BOOST_FUSION_ADAPT_STRUCT(
        smarts_pattern_,
        (atom_, atom)
        (closure_, closure)
        (std::vector<bond_>, bonds)
        (boost::optional<bond_>, last_bond))

/* Convert closure_ind_ into an integer ID */
static
int parse_closure_ind(const closure_ind_& closure_ind) {
    int closure_id;
    if (const char* c = boost::get<char>(&closure_ind)) {
        closure_id = *c - '0';
        if (closure_id < 0 || closure_id > 9)
            MSYS_FAIL("SMARTS BUG: invalid closure id");
    } else if (const bf::vector2<char, char>* v
            = boost::get<bf::vector2<char, char> >(&closure_ind)) {
        if (bf::at_c<0>(*v) < '0' || bf::at_c<0>(*v) > '9'
                || bf::at_c<1>(*v) < '0' || bf::at_c<1>(*v) > '9')
            MSYS_FAIL("SMARTS BUG: invalid closure id");
        closure_id = (bf::at_c<0>(*v)-'0') * 10 + (bf::at_c<1>(*v)-'0');
    } else
        MSYS_FAIL("SMARTS BUG: invalid closure id");
    return closure_id;
}

typedef std::string::const_iterator Iterator;
struct SMARTS_grammar : qi::grammar<Iterator, smarts_pattern_(), ascii::space_type>
{ /* Grammar returns smarts_pattern_ */

    /* Wrap parse_closure_ind function as lazy evaluation function using
     * boost::phoenix */
    struct lazy_closure_ind_impl {
        template<typename A> struct result { typedef int type; };
        int operator()(const closure_ind_& closure) const {
            return parse_closure_ind(closure); }
    };
    boost::phoenix::function<lazy_closure_ind_impl> lazy_closure_ind;

    /* Flags for unsupported SMARTS features */
    bool HAS_DIRECTIONAL_BOND;
    bool HAS_CHIRALITY;
    bool HAS_IMPLICIT_H;
    bool HAS_EXPLICIT_D;
    bool HAS_ISOTOPE;
    bool HAS_ATOMCLASS;
    bool HAS_HYBRIDIZATION;
    void resetFlags() {
        HAS_DIRECTIONAL_BOND = false;
        HAS_CHIRALITY = false;
        HAS_IMPLICIT_H = false;
        HAS_EXPLICIT_D = false;
        HAS_ISOTOPE = false;
        HAS_ATOMCLASS = false;
        HAS_HYBRIDIZATION = false;
    }

    /* Constructor defines SMARTS grammar */
    SMARTS_grammar() : SMARTS_grammar::base_type(smarts_pattern) {

        /* Special case the hydrogen atom; define bonds recursively;
         * bond_expression may be empty; raw_atom, hydrogen_expression,
         * and atom_expression return types must be different to differentiate
         * between them during SMARTS matching */
        smarts_pattern %= 
            ((raw_atom | hydrogen_expression | '[' >> atom_expression >> ']')
                >> closure
                >> *('(' >> bond_expression >> smarts_pattern >> ')')
                >> -(bond_expression >> smarts_pattern));

        /* Interpret double symbols first (where possible) before interpreting
         * single symbols */
        raw_atom %= raw_double_sym | raw_single_sym | qi::char_('R');

        /* For some reason, boost::spirit in boost/1.51.0 does not allow
         * "[H" >> -charge >> "]" here (ambiguous template specification) */
        hydrogen_expression %= "[" >> qi::char_('H') >> -charge >> "]";

        atom_expression %= or_atom_spec % ';';
        or_atom_spec %= and_atom_spec % ',';
        and_atom_spec %= (+not_atom_spec) % '&';
        not_atom_spec %= *qi::char_('!') >> atom_spec;

        /* Ignore isotope and atom class specifications; options with the same
         * return type must be handled by the same if-else branch during
         * SMARTS matching */
        atom_spec %= qi::omit[-isotope] >> (recursive_pattern
            | element
            | atomic_number
            | optional_numeric_property
            | charge
            | implicith
            | explicitD
            | chirality)
        /* To disable support for hybridization, uncomment the below line and
         * remove '^' from character set "XHxvRr^" in definition of
         * optional_numeric_property below */
        //    | hybridization)
            >> qi::omit[-atomclass];

        recursive_pattern %= "$(" >> smarts_pattern >> ")";

        /* Interpet double symbols first (where possible) before interpreting
         * single symbols */
        element %= double_sym | raw_double_sym | single_sym | raw_single_sym;

        atomic_number %= qi::char_("#") >> qi::uint_;
        optional_numeric_property %= qi::char_("XHxvRr^") >> -qi::uint_;
        charge %= qi::char_('+') >> qi::uint_ | qi::char_('-') >> qi::uint_
            | +qi::char_('+') | +qi::char_('-');

        /* Flag implicit h and explicit D; matcher will treat as H and X */
        implicith %= qi::char_("h") >> -qi::uint_
            >> qi::eps[boost::phoenix::ref(HAS_IMPLICIT_H) = true];
        explicitD %= qi::char_("D") >> -qi::uint_
            >> qi::eps[boost::phoenix::ref(HAS_EXPLICIT_D) = true];

        /* Flag isotope and atom class */
        isotope = qi::uint_ >> qi::eps[boost::phoenix::ref(HAS_ISOTOPE) = true];
        atomclass = ':' >> qi::uint_
            >> qi::eps[boost::phoenix::ref(HAS_ATOMCLASS) = true];

        /* Flag chirality and hybridization, replace with '*' atom
         * specification */
        chirality = (qi::char_('@') >> -(
                    qi::char_('@')
                    | qi::lit("TH") >> qi::char_("12")
                    | qi::lit("AL") >> qi::char_("12")
                    | qi::lit("SP") >> qi::char_("123")
                    | qi::lit("TB") >> qi::uint_ 
                    | qi::lit("OH") >> qi::uint_)
                >> -qi::char_('?'))
            [qi::_val = boost::phoenix::val(std::make_pair(0, NONE))]
            >> qi::eps[boost::phoenix::ref(HAS_CHIRALITY) = true];
        hybridization = (qi::char_("^") >> -qi::uint_)
            [qi::_val = boost::phoenix::val(std::make_pair(0, NONE))]
            >> qi::eps[boost::phoenix::ref(HAS_HYBRIDIZATION) = true];

        /* bond_expression may be empty */
        closure %= *(bond_expression >> closure_ind);

        /* Convert to int output */
        closure_ind = (qi::digit | '%' >> qi::digit >> qi::digit)
            [qi::_val = lazy_closure_ind(qi::_1)];

        bond_expression %= -(or_bond_spec % ';');
        or_bond_spec %= and_bond_spec % ',';
        and_bond_spec %= (+bond_spec) % '&';
        bond_spec %= *qi::char_('!')
            >> (qi::char_("=#:@$") | qi::char_('-') | qi::char_('~')
                    | directional_bond);

        /* Flag directional bond, and replace with '-' bond specification */
        directional_bond = (qi::char_("/\\") >> -qi::char_('\?'))
            [qi::_val = boost::phoenix::val('-')]
            >> qi::eps[boost::phoenix::ref(HAS_DIRECTIONAL_BOND) = true];

    }

    /* Define rules in grammar and their return types */
    qi::rule<Iterator, smarts_pattern_(), ascii::space_type>
        smarts_pattern;
    qi::rule<Iterator, raw_atom_(), ascii::space_type> raw_atom;
    qi::rule<Iterator, hydrogen_expression_(), ascii::space_type>
        hydrogen_expression;
    qi::rule<Iterator, atom_expression_(), ascii::space_type>
        atom_expression;
    qi::rule<Iterator, or_atom_spec_(), ascii::space_type> or_atom_spec;
    qi::rule<Iterator, and_atom_spec_(), ascii::space_type> and_atom_spec;
    qi::rule<Iterator, not_atom_spec_(), ascii::space_type> not_atom_spec;
    qi::rule<Iterator, atom_spec_(), ascii::space_type> atom_spec;
    qi::rule<Iterator, smarts_pattern_(), ascii::space_type>
        recursive_pattern;
    qi::rule<Iterator, element_(), ascii::space_type> element;
    qi::rule<Iterator, atomic_number_(), ascii::space_type>
        atomic_number;
    qi::rule<Iterator, optional_numeric_property_(), ascii::space_type>
        optional_numeric_property;
    qi::rule<Iterator, charge_(), ascii::space_type> charge;
    qi::rule<Iterator, optional_numeric_property_(), ascii::space_type>
        implicith;
    qi::rule<Iterator, optional_numeric_property_(), ascii::space_type>
        explicitD;
    qi::rule<Iterator, unsigned(), ascii::space_type> isotope;
    qi::rule<Iterator, unsigned(), ascii::space_type> atomclass;
    qi::rule<Iterator, element_(), ascii::space_type> chirality;
    qi::rule<Iterator, element_(), ascii::space_type> hybridization;
    qi::rule<Iterator, int(), ascii::space_type> closure_ind;
    qi::rule<Iterator, closure_(), ascii::space_type> closure;
    qi::rule<Iterator, bond_expression_(), ascii::space_type>
        bond_expression;
    qi::rule<Iterator, or_bond_spec_(), ascii::space_type> or_bond_spec;
    qi::rule<Iterator, and_bond_spec_(), ascii::space_type> and_bond_spec;
    qi::rule<Iterator, bond_spec_(), ascii::space_type> bond_spec;
    qi::rule<Iterator, char(), ascii::space_type> directional_bond;
};

/******************* Helper functions for SMARTS matching *********************/

static
bool match_element(Id atom, AnnotatedSystemPtr sys, const element_& elem) {
    if (elem.first != 0
            && sys->system()->atom(atom).atomic_number != elem.first)
        return false;
    if (elem.second == NONE)
        return true;
    if (elem.second == AROM && sys->atomAromatic(atom))
        return true;
    if (elem.second == ALIPH && !sys->atomAromatic(atom))
        return true;
    return false;
}

static
bool match_charge(Id atom, AnnotatedSystemPtr sys, const charge_& charge) {
    if (const bf::vector2<char, unsigned>* tmp
            = boost::get<bf::vector2<char, unsigned> >(&charge)) {
        if (bf::at_c<0>(*tmp) == '+')
            return (int(bf::at_c<1>(*tmp))
                    == sys->system()->atom(atom).formal_charge);
        else if (bf::at_c<0>(*tmp) == '-')
            return (-int(bf::at_c<1>(*tmp))
                    == sys->system()->atom(atom).formal_charge);
        else
            MSYS_FAIL("SMARTS BUG; unrecognized charge");
    } else if (const std::vector<char>* tmp
            = boost::get<std::vector<char> >(&charge)) {
        if (tmp->size() == 0)
            MSYS_FAIL("SMARTS BUG; unrecognized charge");
        if (tmp->at(0) == '+')
            /* Assume remaining characters are all '+' */
            return (int(tmp->size())
                    == sys->system()->atom(atom).formal_charge);
        else if (tmp->at(0) == '-')
            /* Assume remaining characters are all '-' */
            return (-int(tmp->size())
                    == sys->system()->atom(atom).formal_charge);
        else
            MSYS_FAIL("SMARTS BUG; unrecognized charge");
    } else
        MSYS_FAIL("SMARTS BUG; unrecognized charge");
}

static
bool match_atom_spec(Id atom, AnnotatedSystemPtr sys, const atom_spec_&
        aspec) {
    if (const SmartsPatternImplPtr* pattern
            = boost::get<SmartsPatternImplPtr>(&aspec)) {
        /* Recursive SMARTS */
        std::vector<IdList> matches;
        return (*pattern)->matchSmartsPattern(sys, atom, matches, true);
    } else if (const element_* elem = boost::get<element_>(&aspec))
        /* Element */
        return match_element(atom, sys, *elem);
    else if (const atomic_number_* anum = boost::get<atomic_number_>(&aspec))
        /* Atomic number */
        return (int(bf::at_c<1>(*anum))
                == sys->system()->atom(atom).atomic_number);
    else if (const charge_* charge = boost::get<charge_>(&aspec))
        /* Charge */
        return match_charge(atom, sys, *charge);
    else if (const optional_numeric_property_* opt
            = boost::get<optional_numeric_property_>(&aspec)) {
        int val = -1;
        if (bf::at_c<1>(*opt))
            val = *(bf::at_c<1>(*opt));
        switch (bf::at_c<0>(*opt)) {
            case 'X':
            case 'D':
                /* Total connections */
                if (val == -1) val = 1;
                return sys->atomDegree(atom) == val;
            case 'H':
            case 'h':
                /* Hydrogen connections */
                if (val == -1) val = 1;
                return sys->atomHcount(atom) == val;
            case 'x':
                /* Ring connections */
                if (val == -1)
                    return sys->atomRingBonds(atom) > 0;
                else
                    return sys->atomRingBonds(atom) == val;
            case 'v':
                /* Total connections counting multiplicity from bond orders */
                if (val == -1) val = 1;
                return sys->atomValence(atom) == val;
            case 'R':
                /* In given number of rings */
                if (val == -1)
                    return sys->atomRingCount(atom);
                else
                    return sys->atomRingCount(atom) == val;
            case 'r':
                /* In ring of given size */
                if (val == -1)
                    return sys->atomRingCount(atom);
                else
                    return sys->atomInRingSize(atom, val);
            case '^':
                if (val == -1) val = 1;
                return sys->atomHybridization(atom) == val;
            default:
                MSYS_FAIL("SMARTS BUG; unrecognized numeric property");
        }
    } else
        /* Includes smarts_pattern_ case, which should have already been
         * converted to SMARTS_pattern type */
        MSYS_FAIL("SMARTS BUG; unrecognized atom spec");
}

static
bool match_not_atom_spec(Id atom, AnnotatedSystemPtr sys,
        const not_atom_spec_& spec) {
    if (bf::at_c<0>(spec).size() % 2 == 0)
        return match_atom_spec(atom, sys, bf::at_c<1>(spec));
    else
        return (!match_atom_spec(atom, sys, bf::at_c<1>(spec)));
}

static
bool match_and_atom_spec(Id atom, AnnotatedSystemPtr sys,
        const and_atom_spec_& spec) {
    for (unsigned i = 0; i < spec.size(); ++i)
        for (unsigned j = 0; j < spec[i].size(); ++j)
            if (!match_not_atom_spec(atom, sys, spec[i][j]))
                return false;
    return true;
}

static
bool match_or_atom_spec(Id atom, AnnotatedSystemPtr sys,
        const or_atom_spec_& spec) {
    for (unsigned i = 0; i < spec.size(); ++i)
        if (match_and_atom_spec(atom, sys, spec[i]))
            return true;
    return false;
}

static
bool match_atom_expression(Id atom, AnnotatedSystemPtr sys,
        const atom_expression_& expr) {
    for (unsigned i = 0; i < expr.size(); ++i)
        if (!match_or_atom_spec(atom, sys, expr[i]))
            return false;
    return true;
}

static
bool match_raw_atom(Id atom, AnnotatedSystemPtr sys,
        const raw_atom_& raw) {
    if (const element_* elem = boost::get<element_>(&raw))
        return match_element(atom, sys, *elem);
    else if (const char* R = boost::get<char>(&raw)) {
        if (*R == 'R')
            return sys->atomRingCount(atom) > 0;
        else
            MSYS_FAIL("SMARTS BUG; unrecognized raw element");
    } else
        MSYS_FAIL("SMARTS BUG; unrecognized raw element");
}

static
bool match_hydrogen_expression(Id atom, AnnotatedSystemPtr sys,
        const hydrogen_expression_& hexpr) {
    if (sys->system()->atom(atom).atomic_number != 1)
        return false;
    if (!bf::at_c<1>(hexpr))
        return true;
    return match_charge(atom, sys, *(bf::at_c<1>(hexpr)));
}

static
bool match_atom(Id atom, AnnotatedSystemPtr sys, const atom_& a) {
    if (const raw_atom_* raw = boost::get<raw_atom_>(&a))
        return match_raw_atom(atom, sys, *raw);
    else if (const hydrogen_expression_* hexpr
            = boost::get<hydrogen_expression_>(&a))
        return match_hydrogen_expression(atom, sys, *hexpr);
    else if (const atom_expression_* expr
            = boost::get<atom_expression_>(&a))
        return match_atom_expression(atom, sys, *expr);
    else
        MSYS_FAIL("SMARTS BUG: unrecognized atom");
}

static
bool match_bond_spec(Id bond, AnnotatedSystemPtr sys,
        const bond_spec_& spec) {
    /* Assume all characters of spec other than the last are '!' */
    bool reverse = (spec.size() % 2 == 1);
    switch (spec[spec.size()-1]) {
        case '~': return (reverse ^ false);
        case '-': return (reverse ^ (sys->system()->bond(bond).order != 1
                              || sys->bondAromatic(bond)));
        case '=': return (reverse ^ (sys->system()->bond(bond).order != 2
                              || sys->bondAromatic(bond)));
        case '#': return (reverse ^ (sys->system()->bond(bond).order != 3
                              || sys->bondAromatic(bond)));
        case '$': return (reverse ^ (sys->system()->bond(bond).order != 4
                              || sys->bondAromatic(bond)));
        case '@': return (reverse ^ (sys->bondRingCount(bond) == 0));
        case ':': return (reverse ^ (!sys->bondAromatic(bond)));
        default: MSYS_FAIL("SMARTS BUG; unrecognized bond spec");
    }
}

static
bool match_and_bond_spec(Id bond, AnnotatedSystemPtr sys,
        const and_bond_spec_& spec) {
    for (unsigned i = 0; i < spec.size(); ++i)
        for (unsigned j = 0; j < spec[i].size(); ++j)
            if (!match_bond_spec(bond, sys, spec[i][j]))
                return false;
    return true;
}

static
bool match_or_bond_spec(Id bond, AnnotatedSystemPtr sys,
        const or_bond_spec_& spec) {
    for (unsigned i = 0; i < spec.size(); ++i)
        if (match_and_bond_spec(bond, sys, spec[i]))
            return true;
    return false;
}

static
bool match_bond_expression(Id bond, AnnotatedSystemPtr sys,
        const bond_expression_& expr) {
    if (expr) {
        for (unsigned i = 0; i < expr->size(); ++i)
            if (!match_or_bond_spec(bond, sys, expr->at(i)))
                return false;
        return true;
    } else if (sys->system()->bond(bond).order == 1
            || sys->bondAromatic(bond))
        /* If bond is unspecified, match single and aromatic bonds */
        return true;
    else
        return false;
}


/****************** Implementation of SmartsPattern class ********************/
SmartsPattern::SmartsPattern(std::string const& pattern)
: _pattern(pattern) {
    std::stringstream ss;
    _impl = SmartsPatternImpl::create(pattern, ss);
    _warnings = ss.str();
}

Id SmartsPattern::atomCount() const {
    return _impl->atomCount();
}

/* Constructor uses grammar to create a recursive smarts_pattern_ object, then
 * converts the smarts_pattern_ into a list of atom and bond expressions in the
 * order they were specified, with closure bonds appearing at the location of
 * their second specification. Matching and checking of closure symbols is done
 * here. */
SmartsPatternImpl::SmartsPatternImpl(std::string const& pattern, 
                                     std::ostream& log) {
    if (pattern.empty()) return;

    /* Use grammar to parse SMARTS pattern */
    SMARTS_grammar grammar;
    grammar.resetFlags();
    std::string::const_iterator iter = pattern.begin();
    std::string::const_iterator end = pattern.end();
    smarts_pattern_ smarts;
    bool success = qi::phrase_parse(iter, end, grammar, ascii::space,
            smarts);
    /* phrase_parse returns success if only part of the string is matched;
     * need to check that iter == end to ensure complete match */
    if (!success || iter != end)
        MSYS_FAIL("Invalid SMARTS string: " + pattern);

    /* Convert to list of atoms and bonds, and check ring closures */
    std::map<int, std::pair<bond_expression_, unsigned> > closure_map;
    if (convertSmarts(smarts, closure_map) && closure_map.size() == 0) {
        /* Warn about unsupported SMARTS features */
        if (grammar.HAS_DIRECTIONAL_BOND)
            log << "...WARNING: replacing directional bond "
                << "specification with \'-\'" << std::endl;
        if (grammar.HAS_CHIRALITY)
            log << "...WARNING: replacing chirality "
                << "specification with wild \'*\'" << std::endl;
        if (grammar.HAS_HYBRIDIZATION)
            log << "...WARNING: replacing hybridization "
                << "specification with wild \'*\'" << std::endl;
        if (grammar.HAS_IMPLICIT_H)
            log << "...WARNING: treating implicit h "
                << "specification as explicit H" << std::endl;
        if (grammar.HAS_EXPLICIT_D)
            log << "...WARNING: treating D "
                << "specification as X" << std::endl;
        if (grammar.HAS_ISOTOPE)
            log << "...WARNING: ignoring isotope "
                << "specification" << std::endl;
        if (grammar.HAS_ATOMCLASS)
            log << "...WARNING: ignoring atom class "
                << "specification" << std::endl;
    } else {
        MSYS_FAIL("Invalid SMARTS string: " + pattern);
    }
}

bool SmartsPatternImpl::convertSmarts(const smarts_pattern_& smarts,
        std::map<int, std::pair<bond_expression_, unsigned> >& closure_map) {
    _atoms.push_back(smarts.atom);
    unsigned aid = _atoms.size() - 1;
    /* Convert smarts_pattern_ into SmartsPattern in recursive SMARTS */
    if (atom_expression_* expr
            = boost::get<atom_expression_>(&_atoms[aid])) {
        for (unsigned i = 0; i < expr->size(); ++i)
        for (unsigned j = 0; j < expr->at(i).size(); ++j)
        for (unsigned k = 0; k < expr->at(i)[j].size(); ++k)
        for (unsigned l = 0; l < expr->at(i)[j][k].size(); ++l) {
            atom_spec_ spec = bf::at_c<1>(expr->at(i)[j][k][l]);
            if (smarts_pattern_* pattern = boost::get<smarts_pattern_>(&spec)) {
                SmartsPatternImplPtr recursive_smarts = create("", std::cout);
                /* Recursive SMARTS operate on a new set of closure indices */
                std::map<int, std::pair<bond_expression_, unsigned> > new_map;
                recursive_smarts->convertSmarts(*pattern, new_map);
                bf::at_c<1>(expr->at(i)[j][k][l]) = recursive_smarts;
            }
        }
    }
    /* Handle closures */
    for (unsigned i = 0; i < smarts.closure.size(); ++i) {
        int ind = bf::at_c<1>(smarts.closure[i]);
        std::map<int, std::pair<bond_expression_, unsigned> >::iterator
            iter = closure_map.find(ind);
        if (iter == closure_map.end()) {
            /* Save this atom and bond expression as a ring opening */
            closure_map.insert(std::make_pair(ind, std::make_pair(
                            bf::at_c<0>(smarts.closure[i]), aid)));
        } else {
            /* Check that closure bond expressions match, add closure bond */
            const bond_expression_& expr = bf::at_c<0>(smarts.closure[i]);
            if (iter->second.first) {
                if (expr && expr != iter->second.first)
                    return false;
                _bonds.push_back(boost::make_tuple(aid,
                            iter->second.first, iter->second.second));
            } else {
                _bonds.push_back(boost::make_tuple(aid, expr,
                            iter->second.second));
            }
            closure_map.erase(iter);
        }
    }
    /* Handle bonds recursively */
    for (unsigned i = 0; i < smarts.bonds.size(); ++i) {
        const bond_expression_& expr = bf::at_c<0>(smarts.bonds[i]);
        _bonds.push_back(boost::make_tuple(aid, expr,
                    _atoms.size()));
        bool success = convertSmarts(bf::at_c<1>(smarts.bonds[i]).get(),
                closure_map);
        if (!success) return false;
    }
    if (smarts.last_bond) {
        const bond_expression_& expr = bf::at_c<0>(*smarts.last_bond);
        _bonds.push_back(boost::make_tuple(aid, expr,
                    _atoms.size()));
        bool success = convertSmarts(bf::at_c<1>(*smarts.last_bond).get(),
                closure_map);
        if (!success) return false;
    }
    return true;
}


MultiIdList SmartsPattern::findMatches(AnnotatedSystemPtr sys,
        IdList const& atoms) const {

    MultiIdList matches;
    BOOST_FOREACH(Id id, atoms) {
        if (sys->system()->atom(id).atomic_number < 1)
            continue;
        _impl->matchSmartsPattern(sys, id, matches);
    }
    return matches;
}

bool SmartsPattern::match(AnnotatedSystemPtr sys) const {
    SystemPtr mol = sys->system();
    MultiIdList matches;
    for (System::iterator i=mol->atomBegin(), e=mol->atomEnd(); i!=e; ++i) {
        if (mol->atomFAST(*i).atomic_number < 1) continue;
        if (_impl->matchSmartsPattern(sys, *i, matches, true)) return true;
    }
    return false;
}

namespace {
    /* Use this with filteredBondsPerAtom to obtain non-pseudo bonds for
     * atoms */
    struct filter_t {
        filter_t(SystemPtr sys) : _sys(sys) { }
        bool operator()(const bond_t& b) const {
            return (_sys->atom(b.i).atomic_number > 0
                    && _sys->atom(b.j).atomic_number > 0);
        }
        SystemPtr _sys;
    };
}

bool SmartsPatternImpl::matchSmartsPattern(AnnotatedSystemPtr sys, Id atom,
        std::vector<IdList>& matches, bool match_single) const {

    if (_atoms.size() == 0)
        return false;
    if (!match_atom(atom, sys, _atoms[0]))
        return false;
    if (_bonds.size() == 0) {
        matches.push_back(IdList(1, atom));
        return true;
    }

    filter_t filter(sys->system());

    if (sys->system()->filteredBondsForAtom(atom, filter).size() == 0)
        return false;

    /* Two-way maps of currently matched atom expressions and atoms in system */
    IdList smarts_to_sys(_atoms.size(), BadId);
    IdList sys_to_smarts(sys->system()->maxAtomId(), BadId);

    /* Stack keeps track of the choice of system bond for each matched bond
     * expression. The next bond expression to match is
     * _bonds[bond_choices.size()-1]. */
    std::stack<unsigned> bond_choices;

    smarts_to_sys[0] = atom;
    sys_to_smarts[atom] = 0;
    bond_choices.push(0);
    bool matched_any = false;

    while (true) {
        const boost::tuple<unsigned, bond_expression_, unsigned>& bond_tuple
            = _bonds[bond_choices.size() - 1];
        Id ai = smarts_to_sys[bond_tuple.get<0>()];
        Id aj = smarts_to_sys[bond_tuple.get<2>()];
        bool closure = (bond_tuple.get<2>() < bond_tuple.get<0>());
        if (ai == BadId)
            MSYS_FAIL("VIPARR_BUG: Bond refers to unmatched atom");
        if (closure && aj == BadId)
            MSYS_FAIL("VIPARR_BUG: Closure bond refers to unmatched atom");
        if (bond_choices.top() >= sys->system()->filteredBondsForAtom(ai,
                    filter).size())
            MSYS_FAIL("SMARTS BUG: Bond choice exceeds number of bonds");
        /* Try matching to the system bond indicated by the top of the stack */
        Id bond = sys->system()->filteredBondsForAtom(ai,
                filter)[bond_choices.top()];
        if (!match_bond_expression(bond, sys, bond_tuple.get<1>())
                || (closure && sys->system()->bond(bond).other(ai) != aj)
                || (!closure
                    && (sys_to_smarts[sys->system()->bond(bond).other(ai)]
                        != BadId || !match_atom(
                            sys->system()->bond(bond).other(ai), sys,
                        _atoms[bond_tuple.get<2>()])))) {
            /* Bond does not match */
            Id top_atom = smarts_to_sys[_bonds[
                bond_choices.size()-1].get<0>()];
            while (bond_choices.top() == sys->system()->filteredBondsForAtom(
                        top_atom, filter).size()-1) {
                /* Have tried all possible system bonds for this atom; pop top
                 * bond choice off of the stack */
                bond_choices.pop();
                if (bond_choices.size() == 0)
                    /* Explored all possible matches from starting atom;
                     * return */
                    return matched_any;
                const boost::tuple<unsigned, bond_expression_, unsigned>&
                    top_tuple = _bonds[bond_choices.size() - 1];
                if (top_tuple.get<2>() > top_tuple.get<0>()) {
                    /* If top bond is not a closure bond, undo the last atom
                     * match */
                    sys_to_smarts[smarts_to_sys[top_tuple.get<2>()]]
                        = BadId;
                    smarts_to_sys[top_tuple.get<2>()] = BadId;
                }
                top_atom = smarts_to_sys[top_tuple.get<0>()];
            }
            /* Try next system bond */
            bond_choices.top() += 1;
        } else {
            /* Bond matches */
            if (!closure) {
                /* If bond is not a closure bond, save the matched atom */
                smarts_to_sys[bond_tuple.get<2>()] = sys->system()->bond(
                        bond).other(ai);
                sys_to_smarts[sys->system()->bond(bond).other(ai)]
                    = bond_tuple.get<2>();
            }
            if (bond_choices.size() == _bonds.size()) {
                /* Found complete match for SMARTS pattern */
                matches.push_back(smarts_to_sys);
                matched_any = true;
                if (match_single)
                    return true;
                if (!closure) {
                    /* If bond is not a closure bond, undo this last match */
                    smarts_to_sys[bond_tuple.get<2>()] = BadId;
                    sys_to_smarts[sys->system()->bond(bond).other(ai)] = BadId;
                }
                Id top_atom = smarts_to_sys[_bonds[
                    bond_choices.size()-1].get<0>()];
                while (bond_choices.top()
                        == sys->system()->filteredBondsForAtom(top_atom,
                            filter).size()-1) {
                    /* Have tried all possible system bonds for this atom; pop
                     * top bond choice off of the stack */
                    bond_choices.pop();
                    if (bond_choices.size() == 0)
                        /* Explored all possible matches from starting atom;
                         * return */
                        return matched_any;
                    const boost::tuple<unsigned, bond_expression_, unsigned>&
                        top_tuple = _bonds[bond_choices.size() - 1];
                    if (top_tuple.get<2>() > top_tuple.get<0>()) {
                        /* If top bond is not a closure bond, undo the last atom
                         * match */
                        sys_to_smarts[smarts_to_sys[top_tuple.get<2>()]]
                            = BadId;
                        smarts_to_sys[top_tuple.get<2>()] = BadId;
                    }
                    top_atom = smarts_to_sys[top_tuple.get<0>()];
                }
                /* Try next system bond */
                bond_choices.top() += 1;
            } else
                /* Not yet a complete match; move on to next bond expression in
                 * SMARTS_pattern */
                bond_choices.push(0);
        }
    }
}
