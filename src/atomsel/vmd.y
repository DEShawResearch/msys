%pure-parser
%name-prefix="vmd_"
%locations
%defines
%error-verbose
%parse-param { VMD* context }
%lex-param { void* scanner }

%{
    #include "vmd.hxx"
    using namespace desres::msys::atomsel;
%}

%union {
    int         ival;
    double      fval;
    char*       sval;
    Regex*      regex;

    Keyword*    key;
    StringPredicate* strfctn;

    IntList*    ilist;
    FloatList*  flist;
    StringList* slist;

    Predicate*  pred;
    Expression* expr;
}

%token <ival> CMP

%token <key> INTKEY FLTKEY STRKEY
%token WITHIN EXWITHIN PBWITHIN WITHINBONDS NEAREST SAME 
%left <key> SINGLE

%token <ival> IVAL
%token <fval> FVAL
%token <sval> SVAL FUNC
%token <regex> REGEX

%token <pred> MACRO 
%token <strfctn> STRFCTN
%nonassoc OF AS

%nonassoc TO
%token ERR

%left OR
%left AND

%left <sval> ADD SUB
%left <sval> MUL DIV MOD
%left <sval> EXP
%left nonassoc NOT

%type <ilist> intlist
%type <flist> fltlist
%type <slist> strlist stringlist

%type <pred>  selection key_value_list
%type <expr>  expression
%type <fval>  number
%type <key>   key

%destructor { free($$); } SVAL
%destructor {delete $$; } intlist fltlist strlist key_value_list 

%{
    int vmd_lex(YYSTYPE* lvalp, YYLTYPE* llocp, void* scanner);
    void vmd_error(YYLTYPE* locp, VMD* context, const char* err) {
      if (context->error.empty()) context->error = err;
    }

    #define scanner context->scanner
%}

%%

start:
      selection { context->result = $1; }
    | error     { context->result = NULL; }
    ;

selection: 
      '(' selection ')'         { $$ = $2;                       }
    | NOT selection             { $$ = context->make_not($2);    }
    | SINGLE                    { $$ = context->make_single($1); }
    | MACRO                     { $$ = $1; }
    | selection AND selection   { $$ = context->make_and($1,$3); }
    | selection OR  selection   { $$ = context->make_or($1,$3);  }
    | key_value_list            { $$ = $1; }
    | WITHIN number OF selection   { $$ = context->make_within($2,$4); }
    | EXWITHIN number OF selection { $$ = context->make_exwithin($2,$4); }
    | PBWITHIN number OF selection { $$ = context->make_pbwithin($2,$4); }
    | WITHINBONDS IVAL OF selection { $$ = context->make_withinbonds($2,$4); }
    | SAME key AS selection     { $$= context->make_same($2,$4); }
    | NEAREST IVAL TO selection { $$ = context->make_nearest($2,$4); }
    | expression CMP expression { $$ = context->make_compare($2,$1,$3); }
    ;

key:
      INTKEY { $$ = $1; }
    | FLTKEY { $$ = $1; }
    | STRKEY { $$ = $1; }
    ;

key_value_list:
      INTKEY intlist            { $$ = context->make_key($1,$2); }
    | FLTKEY fltlist            { $$ = context->make_key($1,$2); }
    | STRKEY strlist            { $$ = context->make_key($1,$2); }
    | STRFCTN stringlist        { $$ = context->make_strfctn($1,$2); }
      ;

number:
      IVAL                      { $$ = $1; }
    | FVAL                      { $$ = $1; }
    ;
 
intlist:
      IVAL                      { $$ = new IntList($1); }
    | IVAL TO IVAL              { $$ = new IntList($1,$3); }
    | intlist IVAL              { $1->add($2); $$=$1; }
    | intlist IVAL TO IVAL      { $1->add($2,$4); $$=$1; }
    ;

fltlist:
      FVAL                      { $$ = new FloatList($1); }
    | FVAL TO FVAL              { $$ = new FloatList($1,$3); }
    | fltlist FVAL              { $1->add($2); $$=$1; }
    | fltlist FVAL TO FVAL      { $1->add($2,$4); $$=$1; }
    ;

strlist:
      SVAL                      { $$ = new StringList($1); free($1); }
    | REGEX                     { $$ = new StringList($1); }
    | IVAL                      { char buf[32];
                                  sprintf(buf, "%d", $1); 
                                  $$ = new StringList(buf); }
    | strlist SVAL              { $1->add($2); free($2); $$=$1;      }
    | strlist REGEX             { $1->add($2); $$=$1;      }
    | strlist IVAL              { char buf[32];
                                  sprintf(buf, "%d", $2); 
                                  $1->add(buf); $$=$1;     }
    ;

stringlist:
      SVAL { $$ = new StringList($1); free($1); }
    | stringlist SVAL { $1->add($2); free($2); $$=$1; }
    ;

expression:
      FVAL                      { $$ = context->make_exp($1); }
    | IVAL                      { $$ = context->make_exp($1); }
    | FLTKEY                    { $$ = context->make_exp($1); }
    | INTKEY                    { $$ = context->make_exp($1); }
    | SUB expression            { $$ = context->make_unaexp("-", $2); }
    | '(' expression ')'        { $$ = $2; }
    | FUNC '(' expression ')'   { $$ = context->make_unaexp($1,$3); free($1); }
    | expression ADD expression { $$ = context->make_binexp("+", $1,$3); }
    | expression SUB expression { $$ = context->make_binexp("-", $1,$3); }
    | expression MUL expression { $$ = context->make_binexp("*", $1,$3); }
    | expression DIV expression { $$ = context->make_binexp("/", $1,$3); }
    | expression MOD expression { $$ = context->make_binexp("%", $1,$3); }
    | expression EXP expression { $$ = context->make_binexp("^", $1,$3); }
    ;

