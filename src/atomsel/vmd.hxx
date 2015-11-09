#ifndef msys_atomsel_context_hxx
#define msys_atomsel_context_hxx

#include "expression.hxx"
#include "msys_keyword.hxx"

namespace desres { namespace msys { namespace atomsel {

    class VMD {
        /* defined in vmd.l */
        void init_scanner();
        void destroy_scanner();
    
        SystemPtr sys;
        const float* pos;
        const double* cell;
        const char* txt = 0;

        /* containers to keep pointers alive */
        std::vector<PredicatePtr> predicates;
        std::vector<ExpressionPtr> expressions;
        std::vector<KeywordPtr> keywords;

        Predicate* add(PredicatePtr p);
        Expression* add(ExpressionPtr e);

    public:
        void* scanner = 0;
        Predicate* result = 0;
        std::string error;
    
        VMD(SystemPtr mol, const float* pos, const double* cell) 
        : sys(mol), pos(pos), cell(cell)
        {}

        inline char getc() { return *txt++; }

        PredicatePtr parse(std::string const& txt);

        /* internal parser interface */
        Keyword* find_key(const char* s);
        Keyword* find_single(const char* s);
        Predicate* find_macro(const char* s);
        StringPredicate* find_strfctn(const char* s);
        Predicate* make_strfctn(StringPredicate* p, StringList* targets);
        bool find_function(const char* s);
        Predicate* make_not(Predicate* sub);
        Predicate* make_and(Predicate* lhs, Predicate* rhs);
        Predicate* make_or(Predicate* lhs, Predicate* rhs);
        Predicate* make_single(Keyword* key);
        Predicate* make_key(Keyword* key, TargetList* targets);
        Predicate* make_within(double r, Predicate* sub);
        Predicate* make_exwithin(double r, Predicate* sub);
        Predicate* make_pbwithin(double r, Predicate* sub);
        Predicate* make_withinbonds(int n, Predicate* sub);
        Predicate* make_same(Keyword* key, Predicate* sub);
        Predicate* make_nearest(int r, Predicate* sub);
        Predicate* make_pbnearest(int r, Predicate* sub);
        Predicate* make_compare(int cmp, Expression* lhs, Expression* rhs);

        Expression* make_exp(Dbl x);
        Expression* make_exp(int x);
        Expression* make_exp(Keyword* key);
        Expression* make_binexp(const char* s, Expression* L, Expression* R);
        Expression* make_unaexp(const char* s, Expression* S);
    };

}}}

#endif
