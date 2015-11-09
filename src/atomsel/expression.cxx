#include "expression.hxx"
#include "keyword.hxx"
#include <functional>
#include <cmath>
#include <stdexcept>
#include <sstream>

using namespace desres::msys::atomsel;

namespace {
    typedef double (*Function)(double);

    double identity(double v) { return v; }
    double negate(double v) { return -v; }
    double sqr(double v) { return v*v; }

    class UnaryExpression : public Expression {
        const std::string sym;
        Function exp;
        ExpressionPtr sub;

        public:
        UnaryExpression(const std::string &_sym, Function _exp,
                ExpressionPtr _sub)
            : sym(_sym), exp(_exp), sub(_sub) {}
        void eval( const Selection& s, std::vector<Dbl>& v ) {
            sub->eval(s,v);
            for (unsigned i=0; i<s.size(); i++) {
                if (s[i]) v[i] = exp(v[i]);
            }
        }
    };

    template <typename T>
        struct Power {
            T operator()(const T& base, const T& exponent) {
                return pow(base,exponent);
            }
        };

    template <class OP>
        class BinaryExpression : public Expression {
            const std::string sym;
            OP op;
            ExpressionPtr lhs;
            ExpressionPtr rhs;

            public:
            BinaryExpression(const std::string &_sym,
                    ExpressionPtr _lhs, ExpressionPtr _rhs )
                : sym(_sym), lhs(_lhs), rhs(_rhs) {}
            void eval( const Selection& s, std::vector<Dbl>& v ) {
                std::vector<Dbl> L(v.size()), R(v.size());
                lhs->eval(s,L);
                rhs->eval(s,R);
                for (unsigned i=0; i<s.size(); i++) {
                    if (s[i]) v[i] = op(L[i], R[i]);
                }
            }
        };

    class LiteralExpression : public Expression {
        const double d;
        public:
        LiteralExpression(double _d) : d(_d) {}
        void eval( const Selection& s, std::vector<Dbl>& v) {
            for (unsigned i=0; i<s.size(); i++) v[i]=d;
        }
    };

    class KeywordExpression : public Expression {
        KeywordPtr key;
        public:
        KeywordExpression(KeywordPtr _key) : key(_key) {}
        void eval( const Selection& s, std::vector<Dbl>& v) {
            switch (key->type) {
                case KEY_INT:
                    {
                        std::vector<Int> vals(s.size());
                        key->iget(s,vals);
                        for (unsigned i=0; i<s.size(); i++) v[i]=vals[i];
                    }
                    break;
                case KEY_DBL:
                    key->dget(s,v);
                    break;
                case KEY_STR:
                    for (unsigned i=0; i<s.size(); i++) v[i]=0;
                    break;
                default:
                    throw std::runtime_error("unknown keyword type");

            }
        }
    };

    template <class OP>
        class RelationPredicate : public Predicate {
            const std::string sym;
            OP op;
            ExpressionPtr lhs, rhs;
            public:
            RelationPredicate( const std::string &_sym,
                    ExpressionPtr _lhs, ExpressionPtr _rhs) 
                : sym(_sym), lhs(_lhs), rhs(_rhs) {}

            void eval( Selection& s ) {
                std::vector<Dbl> L(s.size()), R(s.size());
                lhs->eval(s,L);
                rhs->eval(s,R);
                for (unsigned i=0; i<s.size(); i++) {
                    s[i] &= op(L[i],R[i]);
                }
            }
        };
}

namespace desres { namespace msys { namespace atomsel {

    ExpressionPtr unary_expression( const std::string& op, ExpressionPtr sub ) {
        Function exp = 
            op=="+" ? identity :
            op=="-" ? negate :
            op=="sqr" ? sqr : 
            op=="sqrt" ? sqrt :
            op=="abs" ? fabs :
            NULL;
        if (!exp) {
            std::stringstream ss;
            ss << "Unsupported unary operand '" << op << "'";
            throw std::runtime_error(ss.str());
        }
        return ExpressionPtr(new UnaryExpression(op,exp,sub));
    }

    ExpressionPtr binary_expression( const std::string& op, 
            ExpressionPtr lhs, ExpressionPtr rhs ) {
        ExpressionPtr p;
        if (op=="+") 
            p.reset(new BinaryExpression<std::plus<Dbl> >(op,lhs,rhs));
        else if (op=="-")
            p.reset(new BinaryExpression<std::minus<Dbl> >(op,lhs,rhs));
        else if (op=="*")
            p.reset(new BinaryExpression<std::multiplies<Dbl> >(op,lhs,rhs));
        else if (op=="/")
            p.reset(new BinaryExpression<std::divides<Dbl> >(op,lhs,rhs));
        else if (op=="%")
            p.reset(new BinaryExpression<std::modulus<Int> >(op,lhs,rhs));
        else if (op=="**" || op=="^")
            p.reset(new BinaryExpression<Power<Dbl> >(op,lhs,rhs));
        else {
            std::stringstream ss;
            ss << "Unsupported binary operand '" << op << "'";
            throw std::runtime_error(ss.str());
        }
        return p;
    }

    ExpressionPtr literal_expression( double v ) {
        return ExpressionPtr(new LiteralExpression(v));
    }

    ExpressionPtr keyword_expression( KeywordPtr key ) {
        return ExpressionPtr(new KeywordExpression(key));
    }

    PredicatePtr relation_predicate( RelationType op,
            ExpressionPtr lhs, ExpressionPtr rhs ) {
        PredicatePtr p;
        if      (op==RELATION_LT)
            p.reset(new RelationPredicate<std::less<Dbl> >("<",lhs,rhs));
        else if (op==RELATION_GT)
            p.reset(new RelationPredicate<std::greater<Dbl> >(">",lhs,rhs));
        else if (op==RELATION_LE)
            p.reset(new RelationPredicate<std::less_equal<Dbl> >("<=",lhs,rhs));
        else if (op==RELATION_GE)
            p.reset(new RelationPredicate<std::greater_equal<Dbl> >(">=",lhs,rhs));
        else if (op==RELATION_EQ)
            p.reset(new RelationPredicate<std::equal_to<Dbl> >("==",lhs,rhs));
        else if (op==RELATION_NE)
            p.reset(new RelationPredicate<std::not_equal_to<Dbl> >("!=",lhs,rhs));
        else {
            std::stringstream ss;
            ss << "Unsupported relation operand '" << op << "'";
            throw std::runtime_error(ss.str());
        }
        return p;
    }

}}}

