#ifndef msys_atomsel_tokenizer
#define msys_atomsel_tokenizer

#include <stdio.h>
#include <string>
#include "atomsel.h"
#include "selection.hxx"
#include "../system.hxx"

namespace desres { namespace msys { namespace atomsel {

struct Predicate;
struct Query;
struct Token {
    enum RelOp {
        EQ, NE, LT, LE, GE, GT
    };

    int ival;
    double fval;
    const char* sval;
    int slen;
    double (*func)(double);

    std::string str() const { return std::string(sval,sval+slen); }
};

class Tokenizer {
    const char* const s;
    int loc;

public:
    Tokenizer(const char* sel)
    : s(sel), loc(0) {}

    int next(Token* t, System* mol);
    int location() const { return loc; }
};

struct Valist {
    std::vector<int> ival;
    std::vector<double> fval;
    std::vector<std::string> sval;

    std::vector<std::pair<int,int>> irng;
    std::vector<std::pair<double,double>> frng;
    std::vector<std::string> regex;

    void add(int i) { ival.push_back(i); }
    void add(double f) { fval.push_back(f); }
    void add(std::string&& s) { sval.emplace_back(s); }

    void add(int i1, int i2) { irng.emplace_back(i1,i2); }
    void add(double f1, double f2) { frng.emplace_back(f1,f2); }
    void add_regex(std::string&& s) { regex.emplace_back(s); }
};

struct Predicate {
    virtual ~Predicate() = default;
    virtual void eval(Selection& s) = 0;
};

struct BoolPredicate : Predicate {
    System* mol;
    std::string name;
    
    BoolPredicate(System* m, std::string&& s)
    : mol(m), name(s) {}
    virtual void eval(Selection& s);
};

struct KeyPredicate : Predicate {
    Query* q;
    std::string name;
    std::unique_ptr<Valist> va;

    KeyPredicate(Query* q, std::string&& s, Valist* v)
    : q(q), name(s), va(v) {}
    virtual void eval(Selection& s);
};

struct AndPredicate : Predicate {
    std::unique_ptr<Predicate> lhs, rhs;
    AndPredicate(Predicate* L, Predicate* R) : lhs(L), rhs(R) {}
    virtual void eval(Selection& s) {
        lhs->eval(s);
        rhs->eval(s);
    }
};
struct OrPredicate : Predicate {
    std::unique_ptr<Predicate> lhs, rhs;
    OrPredicate(Predicate* L, Predicate* R) : lhs(L), rhs(R) {}
    virtual void eval(Selection& s) {
        Selection s2(s);
        lhs->eval(s);
        rhs->eval(s2);
        s.add(s2);
    }
};
struct NotPredicate : Predicate {
    std::unique_ptr<Predicate> sub;
    NotPredicate(Predicate* S) : sub(S) {}
    virtual void eval(Selection& s) {
        Selection s2(s);
        sub->eval(s2);
        s.subtract(s2);
    }
};
class WithinPredicate : public Predicate {
    System* sys;
    const float* pos;
    const double* cell;
    const float rad;
    std::unique_ptr<Predicate> sub;
    const bool exclude;
    const bool periodic;

public:
    WithinPredicate( System* e, const float* pos, const double* cell, float r, bool excl, bool per, Predicate* s )
    : sys(e), pos(pos), cell(cell), rad(r), sub(s), exclude(excl), periodic(per) {}

  void eval( Selection& s );
};

class WithinBondsPredicate : public Predicate {
  System* sys;
  const int N;
  std::unique_ptr<Predicate> sub;

public:
  WithinBondsPredicate( System* e, int n, Predicate* s )
    : sys(e), N(n), sub(s) {}

  void eval( Selection& s );
};
class KNearestPredicate : public Predicate {
  System* _sys;
  const float* pos;
  const double* cell;
  const unsigned _N;
  const bool periodic;
  std::unique_ptr<Predicate> _sub;

public:
  KNearestPredicate(System* sys, const float* pos, const double* cell, unsigned k, bool per, Predicate* sub)
  : _sys(sys), pos(pos), cell(cell), _N(k), periodic(per), _sub(sub) {}

  void eval(Selection& s);
};

struct SamePredicate : Predicate {
    Query* q;
    std::string name;
    std::unique_ptr<Predicate> sub;

    SamePredicate(Query* q, std::string&& name, Predicate* p)
    : q(q), name(name), sub(p) {}

    void eval( Selection& s );
};

struct Expression {
    virtual ~Expression() = default;
    virtual void eval(Selection const& s, std::vector<double>& v) = 0;
};

struct LitExpr : Expression {
    const double d;
    LitExpr(double d) : d(d) {}
    void eval(Selection const& s, std::vector<double>& v) {
        std::fill(v.begin(), v.end(), d);
    }
};

struct KeyExpr : Expression {
    Query* q;
    std::string name;
    KeyExpr(Query* q, std::string&& s) : q(q), name(s) {}
    void eval(Selection const& s, std::vector<double>& v);
};

struct FuncExpr : Expression {
    double (*func)(double);
    std::unique_ptr<Expression> sub;
    FuncExpr(double (*f)(double), Expression* e)
    : func(f), sub(e) {}
    void eval(Selection const& s, std::vector<double>& v);
};

struct NegExpr : Expression {
    std::unique_ptr<Expression> sub;
    NegExpr(Expression* e) : sub(e) {}
    void eval(Selection const& s, std::vector<double>& v);
};

struct BinExpr : Expression {
    int op;
    std::unique_ptr<Expression> lhs;
    std::unique_ptr<Expression> rhs;
    
    BinExpr(int op, Expression* L, Expression* R)
    : op(op), lhs(L), rhs(R) {}
    void eval(Selection const& s, std::vector<double>& v);
};

struct CmpPredicate : Predicate {
    int cmp;
    std::unique_ptr<Expression> lhs, rhs;

    CmpPredicate(int c, Expression* L, Expression* R) 
    : cmp(c), lhs(L), rhs(R) {}

    void eval( Selection& s );
};

struct Query {
    System* mol = nullptr;
    const float* pos = nullptr;
    const double* cell = nullptr;
    Id id = BadId;  // hack for key selections
    std::unique_ptr<Predicate> pred;

    void parse(std::string const& selection);
};

bool is_keyword(std::string const& name, System* mol);


Selection full_selection(System* sys);

}}}


#endif

