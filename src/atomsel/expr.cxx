#include "token.hxx"
#include <unordered_map>

using namespace desres::msys;
using namespace desres::msys::atomsel;

void CmpPredicate::eval(Selection& s) {
    std::vector<double> L(s.size()), R(s.size());
    lhs->eval(s,L);
    rhs->eval(s,R);
    switch (cmp) {
    case Token::EQ:
        for (Id i=0, n=s.size(); i<n; i++) if (s[i]) s[i]=L[i]==R[i]; break;
    case Token::LT:
        for (Id i=0, n=s.size(); i<n; i++) {
            if (s[i]) s[i]=L[i]<R[i];
        }
        break;
    case Token::GT:
        for (Id i=0, n=s.size(); i<n; i++) if (s[i]) s[i]=L[i]>R[i]; break;
    case Token::LE:
        for (Id i=0, n=s.size(); i<n; i++) if (s[i]) s[i]=L[i]<=R[i]; break;
    case Token::GE:
        for (Id i=0, n=s.size(); i<n; i++) if (s[i]) s[i]=L[i]>=R[i]; break;
    default:;
    }
}

void NegExpr::eval(Selection const& s, std::vector<double>& v) {
    sub->eval(s,v);
    for (auto& x : v) x=-x;
}

void FuncExpr::eval(Selection const& s, std::vector<double>& v) {
    sub->eval(s,v);
    for (Id i=0, n=s.size(); i<n; i++) {
        if (s[i]) v[i] = func(v[i]);
    }
}

void BinExpr::eval(Selection const& s, std::vector<double>& v) {
    std::vector<double> L(v.size()), R(v.size());
    lhs->eval(s,L);
    rhs->eval(s,R);
    switch (op) {
    case ADD: for (Id i=0,n=s.size(); i<n; i++) if (s[i]) v[i]=L[i]+R[i]; break;
    case SUB: for (Id i=0,n=s.size(); i<n; i++) if (s[i]) v[i]=L[i]-R[i]; break;
    case MUL: for (Id i=0,n=s.size(); i<n; i++) if (s[i]) v[i]=L[i]*R[i]; break;
    case DIV: for (Id i=0,n=s.size(); i<n; i++) if (s[i]) v[i]=L[i]/R[i]; break;
    case MOD: for (Id i=0,n=s.size(); i<n; i++) {
                  if (s[i]) v[i]=int(L[i]) % int(R[i]);
              }
              break;
    case EXP: for (Id i=0,n=s.size(); i<n; i++) if (s[i]) v[i]=pow(L[i],R[i]); break;
    }
}

