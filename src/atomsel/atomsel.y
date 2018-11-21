
%name atomselParse

%include {
#include <assert.h>
#include "token.hxx"
using namespace desres::msys::atomsel;

#define YYNOERRORRECOVERY
}

%extra_argument {Query* query}
%token_type {Token}

%nonassoc OF AS.
%nonassoc TO.
%left OR.
%left AND.
%left ADD SUB.
%left MUL DIV MOD.
%left EXP.
%right NEG.
%left NOT.

%type expr {Expression*}
%type pexpr {Expression*}
%type nexpr {Expression*}
%type binexpr {Expression*}
%type selection {Predicate*}
%type num double
%type integer int
%type flt double
%type list    {Valist*}
%type intlist {Valist*}
%type fltlist {Valist*}
%type strlist {Valist*}

%destructor expr { (void)query; delete $$; }
%destructor nexpr { (void)query; delete $$; }
%destructor pexpr { (void)query; delete $$; }
%destructor selection { delete $$; }
%destructor list { delete $$; }
%destructor intlist { delete $$; }
%destructor fltlist { delete $$; }
%destructor strlist { delete $$; }

%syntax_error { query->mol = nullptr; }

input ::= selection(s). { query->pred.reset(s); }
input ::= .

    //WithinPredicate( System* e, const float* pos, const double* cell, float r, bool excl, bool per, Predicate* s )

selection(S) ::= VAL(V).          { S=new BoolPredicate(query->mol,V.str());  }
selection(S) ::= KEY(V) list(v).  { S=new KeyPredicate(query,V.str(),v); }
selection(S) ::= selection(a) AND selection(b). {S=new AndPredicate(a,b); }
selection(S) ::= selection(a) OR selection(b).  {S=new OrPredicate(a,b); }
selection(S) ::= LPAREN selection(s) RPAREN.    {S=s; }
selection(S) ::= MACRO.                         {S=query->pred.release(); }
selection(S) ::= NOT selection(s).              {S=new NotPredicate(s); }
selection(S) ::= WITHIN num(n) OF selection(s).   {S=new WithinPredicate(query->mol, query->pos, NULL, n, false, false, s); }
selection(S) ::= EXWITHIN num(n) OF selection(s). {S=new WithinPredicate(query->mol, query->pos, NULL, n,  true, false, s); }
selection(S) ::= PBWITHIN num(n) OF selection(s). {S=new WithinPredicate(query->mol, query->pos, query->cell, n, false,  true, s); }
selection(S) ::= NEAREST INT(v) TO selection(s).   {S=new KNearestPredicate(query->mol,query->pos, NULL, v.ival, false, s); }
selection(S) ::= WITHINBONDS INT(v) OF selection(s).   {S=new WithinBondsPredicate(query->mol,v.ival, s); }
selection(S) ::= PBNEAREST INT(v) TO selection(s). {S=new KNearestPredicate(query->mol,query->pos,query->cell, v.ival,  true, s); }
selection(S) ::= SAME KEY(v) AS selection(s).   {S=new SamePredicate(query,v.str(),s); }
selection(S) ::= expr(a) CMP(c) expr(b).      {S=new CmpPredicate(c.ival,a,b);}

list(V) ::= integer(T). { V=new Valist; V->add(T); }
list(V) ::= integer(T1) TO integer(T2). { V=new Valist; V->add(T1,T2); }
list(V) ::= VAL(T).     { V=new Valist; V->add(T.str()); }
list(V) ::= REGEX(T).   { V=new Valist; V->add_regex(T.str()); }
list(V) ::= flt(T).     { V=new Valist; V->add(T); }
list(V) ::= flt(T1) TO flt(T2). { V=new Valist; V->add(T1,T2); }

list(V) ::= list(v) integer(T). { V=v; V->add(T); }
list(V) ::= list(v) integer(T1) TO integer(T2). { V=v; V->add(T1,T2); }
list(V) ::= list(v) VAL(T).     { V=v; V->add(T.str()); }
list(V) ::= list(v) REGEX(T).   { V=v; V->add_regex(T.str()); }
list(V) ::= list(v) flt(T).     { V=v; V->add(T); }
list(V) ::= list(v) flt(T1) TO flt(T2). { V=v; V->add(T1,T2); }

integer(I) ::= INT(T).     { I= T.ival; }
integer(I) ::= NEG INT(T). { I=-T.ival; }

flt(I) ::= FLT(T).       { I= T.fval; }
flt(I) ::= NEG FLT(T).   { I=-T.fval; }

num(N) ::= INT(T).      { N=T.ival;  }
num(N) ::= FLT(T).      { N=T.fval;  }
num(N) ::= NEG INT(T).  { N=-T.ival; }
num(N) ::= NEG FLT(T).  { N=-T.fval; }
num(N) ::= PLUS INT(T).  { N=T.ival; }
num(N) ::= PLUS FLT(T).  { N=T.fval; }

expr(E) ::= num(v).         { E=new LitExpr(v); }
expr(E) ::= nexpr(e).       { E=e; }
expr(E) ::= NEG nexpr(e).   { E=new NegExpr(e); }

nexpr(E) ::= KEY(v).                        { E=new KeyExpr(query,v.str()); }
nexpr(E) ::= LPAREN expr(e) RPAREN.         { E=e; }
nexpr(E) ::= FUNC(v) LPAREN expr(e) RPAREN. { E=new FuncExpr(v.func, e); }

expr(E) ::= expr(a) ADD expr(b). {E=new BinExpr(ADD,a,b); }
expr(E) ::= expr(a) SUB expr(b). {E=new BinExpr(SUB,a,b); }
expr(E) ::= expr(a) MUL expr(b). {E=new BinExpr(MUL,a,b); }
expr(E) ::= expr(a) DIV expr(b). {E=new BinExpr(DIV,a,b); }
expr(E) ::= expr(a) MOD expr(b). {E=new BinExpr(MOD,a,b); }
expr(E) ::= expr(a) EXP expr(b). {E=new BinExpr(EXP,a,b); }

// vim: filetype=c
