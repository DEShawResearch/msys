%pure-parser
%name-prefix="msys_smiles_"
%locations
%defines
%error-verbose
%parse-param { Smiles* context }
%lex-param { void* scanner }

%{
#include "smiles.hxx"
#include "../elements.hxx"
using namespace desres::msys::smiles;
using desres::msys::ElementForAbbreviation;
%}

%union {
    atom_t* atom;
    chain_t* chain;
    branch_t *branch;
    ringbond_t *ringbond;
    char    sym[4];
    char    bond;
    int     num;
};

%token<num> ISOTOPE HCOUNT CHARGE CLASS CHIRAL
%token<sym> ATOM SYMBOL
%token<bond> SBOND
%token<ringbond> RINGBOND
%token ERR LPAREN RPAREN DOT PERCENT LBRACKET RBRACKET

%type<atom> atom branched_atom 
%type<atom> bracket_atom isotope symbol chiral hcount charge class rbracket
%type<chain> chain
%type<branch> branch branchlist branchspec
%type<ringbond> ringbondlist 

%{
    int msys_smiles_lex(YYSTYPE* lvalp, YYLTYPE* llocp, void* scanner);
    void msys_smiles_error(YYLTYPE* locp, Smiles* context, const char* err) {
        if (context->error.empty()) context->error = err;
    }
    #define scanner context->scanner
%}

%%

start: 
        chain { context->finish($1); }
      ;
        
chain:
        branched_atom { 
            $$ = context->makeChain($1);
        }
      | chain branched_atom { 
            context->addBond($2->next=$1->last,$2,'-'); 
            $1->last = $2;
            $$ = $1; 
        }
      | chain SBOND branched_atom {
            context->addBond($3->next=$1->last,$3,$2); 
            $1->last = $3;
            $$ = $1; 
        }
      | chain DOT  branched_atom { 
            context->addBond($3->next=$1->last,$3,'.');
            $1->last = $3;
            $$ = $1; 
        }
      ;

branched_atom:
        atom                { $$ = $1; }
      | atom ringbondlist   { 
            $$ = $1; 
            context->addRing($1,$2);
        }
      | atom branchlist { 
            $$ = $1; 
            context->addBranch($1,$2);
        }
      | atom branchlist ringbondlist { 
            $$ = $1;
            context->addBranch($1,$2);
            context->addRing($1,$3);
        }
      | atom ringbondlist branchlist { 
            $$ = $1;
            context->addRing($1,$2);
            context->addBranch($1,$3);
        }
      ;

branchlist:
        branch              { $$ = $1; }
      | branchlist branch   { $2->next = $1; $$ = $2; }
      ;

branch:
        LPAREN branchspec RPAREN { $$ = $2; }
      ; 

branchspec:
        chain               { $$ = context->makeBranch('-', $1); }
      | SBOND chain         { $$ = context->makeBranch($1 , $2); }
      | DOT chain           { $$ = context->makeBranch('.', $2);}
      ;

ringbondlist:
        RINGBOND                { $$ = $1;                  }
      | ringbondlist RINGBOND   { $2->next = $1; $$ = $2;   }
      ;

atom:
        ATOM    { $$ = context->makeAtom(); $$->setName($1); context->addAtom($$,true); }
      | bracket_atom    { $$ = $1;                  context->addAtom($$,false);}
      ;

bracket_atom:
        LBRACKET isotope    { $$ = $2; }
      ;

isotope:
        symbol          { $$ = $1;                  }
      | ISOTOPE symbol  { $$ = $2; $$->isotope = $1;}
      ;

symbol: 
        SYMBOL chiral   { $$ = $2; $$->setName($1); }
      ;

chiral:
        hcount          { $$ = $1;                  }
      | CHIRAL hcount   { $$ = $2; $$->chiral = $1; }
      ;

hcount:
        charge          { $$ = $1;                  }
      | HCOUNT charge   { $$ = $2; $$->hcount = $1; }
      ;

charge:
        class           { $$ = $1;                  }
      | CHARGE class    { $$ = $2; $$->charge = $1; }
      ;

class:
        rbracket        { $$ = $1;                  }
      | CLASS rbracket  { $$ = $2; $$->klass = $1;  }
      ;

rbracket:
        RBRACKET        { $$ = context->makeAtom(); }
      ;

