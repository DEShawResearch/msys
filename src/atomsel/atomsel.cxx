/* Driver template for the LEMON parser generator.
** The author disclaims copyright to this source code.
*/
/* First off, code is included that follows the "include" declaration
** in the input grammar file. */
#include <stdio.h>
#line 4 "atomsel.y"

#include <assert.h>
#include "token.hxx"
using namespace desres::msys::atomsel;

#define YYNOERRORRECOVERY
#line 15 "atomsel.c"
/* Next is all token values, in a form suitable for use by makeheaders.
** This section will be null unless lemon is run with the -m switch.
*/
/* 
** These constants (all generated automatically by the parser generator)
** specify the various kinds of tokens (terminals) that the parser
** understands. 
**
** Each symbol here is a terminal symbol in the grammar.
*/
/* Make sure the INTERFACE macro is defined.
*/
#ifndef INTERFACE
# define INTERFACE 1
#endif
/* The next thing included is series of defines which control
** various aspects of the generated parser.
**    YYCODETYPE         is the data type used for storing terminal
**                       and nonterminal numbers.  "unsigned char" is
**                       used if there are fewer than 250 terminals
**                       and nonterminals.  "int" is used otherwise.
**    YYNOCODE           is a number of type YYCODETYPE which corresponds
**                       to no legal terminal or nonterminal number.  This
**                       number is used to fill in empty slots of the hash 
**                       table.
**    YYFALLBACK         If defined, this indicates that one or more tokens
**                       have fall-back values which should be used if the
**                       original value of the token will not parse.
**    YYACTIONTYPE       is the data type used for storing terminal
**                       and nonterminal numbers.  "unsigned char" is
**                       used if there are fewer than 250 rules and
**                       states combined.  "int" is used otherwise.
**    atomselParseTOKENTYPE     is the data type used for minor tokens given 
**                       directly to the parser from the tokenizer.
**    YYMINORTYPE        is the data type used for all minor tokens.
**                       This is typically a union of many types, one of
**                       which is atomselParseTOKENTYPE.  The entry in the union
**                       for base tokens is called "yy0".
**    YYSTACKDEPTH       is the maximum depth of the parser's stack.  If
**                       zero the stack is dynamically sized using realloc()
**    atomselParseARG_SDECL     A static variable declaration for the %extra_argument
**    atomselParseARG_PDECL     A parameter declaration for the %extra_argument
**    atomselParseARG_STORE     Code to store %extra_argument into yypParser
**    atomselParseARG_FETCH     Code to extract %extra_argument from yypParser
**    YYNSTATE           the combined number of states.
**    YYNRULE            the number of rules in the grammar
**    YYERRORSYMBOL      is the code number of the error symbol.  If not
**                       defined, then do no error processing.
*/
#define YYCODETYPE unsigned char
#define YYNOCODE 47
#define YYACTIONTYPE unsigned char
#define atomselParseTOKENTYPE Token
typedef union {
  int yyinit;
  atomselParseTOKENTYPE yy0;
  Valist* yy14;
  Predicate* yy16;
  Expression* yy19;
  int yy52;
  double yy76;
} YYMINORTYPE;
#ifndef YYSTACKDEPTH
#define YYSTACKDEPTH 100
#endif
#define atomselParseARG_SDECL Query* query;
#define atomselParseARG_PDECL ,Query* query
#define atomselParseARG_FETCH Query* query = yypParser->query
#define atomselParseARG_STORE yypParser->query = query
#define YYNSTATE 102
#define YYNRULE 51
#define YY_NO_ACTION      (YYNSTATE+YYNRULE+2)
#define YY_ACCEPT_ACTION  (YYNSTATE+YYNRULE+1)
#define YY_ERROR_ACTION   (YYNSTATE+YYNRULE)

/* The yyzerominor constant is used to initialize instances of
** YYMINORTYPE objects to zero. */
static const YYMINORTYPE yyzerominor = { 0 };

/* Define the yytestcase() macro to be a no-op if is not already defined
** otherwise.
**
** Applications can choose to define yytestcase() in the %include section
** to a macro that can assist in verifying code coverage.  For production
** code the yytestcase() macro should be turned off.  But it is useful
** for testing.
*/
#ifndef yytestcase
# define yytestcase(X)
#endif


/* Next are the tables used to determine what action to take based on the
** current state and lookahead token.  These tables are used to implement
** functions that take a state number and lookahead value and return an
** action integer.  
**
** Suppose the action integer is N.  Then the action is determined as
** follows
**
**   0 <= N < YYNSTATE                  Shift N.  That is, push the lookahead
**                                      token onto the stack and goto state N.
**
**   YYNSTATE <= N < YYNSTATE+YYNRULE   Reduce by rule N-YYNSTATE.
**
**   N == YYNSTATE+YYNRULE              A syntax error has occurred.
**
**   N == YYNSTATE+YYNRULE+1            The parser accepts its input.
**
**   N == YYNSTATE+YYNRULE+2            No such action.  Denotes unused
**                                      slots in the yy_action[] table.
**
** The action table is constructed as a single large table named yy_action[].
** Given state S and lookahead X, the action is computed as
**
**      yy_action[ yy_shift_ofst[S] + X ]
**
** If the index value yy_shift_ofst[S]+X is out of range or if the value
** yy_lookahead[yy_shift_ofst[S]+X] is not equal to X or if yy_shift_ofst[S]
** is equal to YY_SHIFT_USE_DFLT, it means that the action is not in the table
** and that yy_default[S] should be used instead.  
**
** The formula above is for computing the action when the lookahead is
** a terminal symbol.  If the lookahead is a non-terminal (as occurs after
** a reduce action) then the yy_reduce_ofst[] array is used in place of
** the yy_shift_ofst[] array and YY_REDUCE_USE_DFLT is used in place of
** YY_SHIFT_USE_DFLT.
**
** The following are the tables generated in this section:
**
**  yy_action[]        A single table containing all actions.
**  yy_lookahead[]     A table containing the lookahead for each entry in
**                     yy_action.  Used to detect hash collisions.
**  yy_shift_ofst[]    For each state, the offset into yy_action for
**                     shifting terminals.
**  yy_reduce_ofst[]   For each state, the offset into yy_action for
**                     shifting non-terminals after a reduce.
**  yy_default[]       Default action for each state.
*/
static const YYACTIONTYPE yy_action[] = {
 /*     0 */   103,   13,   15,   17,   18,   19,   20,   17,   18,   19,
 /*    10 */    20,   93,   23,    4,   72,   21,    3,   57,   87,   24,
 /*    20 */    25,   26,   60,   89,   62,   64,   66,   50,   52,   90,
 /*    30 */    48,   68,   23,    4,   72,   21,    3,   73,   87,   24,
 /*    40 */    25,   26,   60,   89,   62,   64,   66,   51,   49,   90,
 /*    50 */    48,   68,   13,   15,   17,   18,   19,   20,   74,   89,
 /*    60 */    54,   55,   22,   99,   23,   90,   48,   94,   14,   40,
 /*    70 */    53,   82,   32,   12,   98,   89,   38,   97,    2,    1,
 /*    80 */    74,   90,   48,   68,  154,   83,   79,   79,   13,   15,
 /*    90 */    17,   18,   19,   20,   58,   13,   15,   17,   18,   19,
 /*   100 */    20,   59,   78,   32,   81,   98,   99,   85,   97,   12,
 /*   110 */    40,   84,   76,   13,   15,   17,   18,   19,   20,   27,
 /*   120 */    35,   74,   98,   75,  100,   97,   77,   79,   94,   14,
 /*   130 */    32,   80,   98,  102,   56,   97,   91,    2,    1,    1,
 /*   140 */    28,   31,   92,   98,   68,   39,   97,   32,   37,   98,
 /*   150 */    98,   88,   97,   97,   29,   32,   30,   98,    5,   41,
 /*   160 */    97,   32,   33,   98,   98,   42,   97,   97,   32,   36,
 /*   170 */    98,   98,   43,   97,   97,   32,    6,   98,    7,   44,
 /*   180 */    97,   32,   34,   98,   98,   45,   97,   97,   61,   32,
 /*   190 */    63,   98,    8,   46,   97,   32,   69,   98,   98,   47,
 /*   200 */    97,   97,   70,   71,   98,   98,    9,   97,   97,  101,
 /*   210 */    75,   98,    2,    1,   97,   95,   80,   91,   65,   10,
 /*   220 */    11,   96,   67,   92,   20,   86,   16,
};
static const YYCODETYPE yy_lookahead[] = {
 /*     0 */     0,    6,    7,    8,    9,   10,   11,    8,    9,   10,
 /*    10 */    11,   35,   12,   13,   14,   15,   16,   38,   18,   19,
 /*    20 */    20,   21,   22,   23,   24,   25,   26,   39,   40,   29,
 /*    30 */    30,   31,   12,   13,   14,   15,   16,   39,   18,   19,
 /*    40 */    20,   21,   22,   23,   24,   25,   26,   12,   12,   29,
 /*    50 */    30,   31,    6,    7,    8,    9,   10,   11,   23,   23,
 /*    60 */    39,   40,   41,   17,   12,   29,   30,   15,   16,   12,
 /*    70 */    12,   14,   33,   27,   35,   23,   37,   38,    4,    5,
 /*    80 */    23,   29,   30,   31,   45,   28,   29,   29,    6,    7,
 /*    90 */     8,    9,   10,   11,   38,    6,    7,    8,    9,   10,
 /*   100 */    11,   38,   40,   33,   39,   35,   17,   37,   38,   27,
 /*   110 */    12,   40,   14,    6,    7,    8,    9,   10,   11,    3,
 /*   120 */    33,   23,   35,   23,   17,   38,   28,   29,   15,   16,
 /*   130 */    33,   29,   35,    0,   37,   38,   23,    4,    5,    5,
 /*   140 */     3,   33,   29,   35,   31,   37,   38,   33,   33,   35,
 /*   150 */    35,   37,   38,   38,    3,   33,    3,   35,    1,   37,
 /*   160 */    38,   33,   33,   35,   35,   37,   38,   38,   33,   33,
 /*   170 */    35,   35,   37,   38,   38,   33,    1,   35,    1,   37,
 /*   180 */    38,   33,   33,   35,   35,   37,   38,   38,   23,   33,
 /*   190 */    23,   35,    3,   37,   38,   33,   33,   35,   35,   37,
 /*   200 */    38,   38,   33,   33,   35,   35,    1,   38,   38,   33,
 /*   210 */    23,   35,    4,    5,   38,   23,   29,   23,   23,    3,
 /*   220 */     2,   29,   15,   29,   11,   17,   16,
};
#define YY_SHIFT_USE_DFLT (-6)
#define YY_SHIFT_MAX 71
static const short yy_shift_ofst[] = {
 /*     0 */     0,   20,   20,   20,   20,   20,   20,   20,   20,   20,
 /*    10 */    20,   20,   52,   52,   52,   52,   52,   52,   52,   52,
 /*    20 */    52,   57,   98,  113,   36,   36,   36,   35,   58,   35,
 /*    30 */    58,   46,   82,   89,  107,   -5,   -1,   -1,  133,  208,
 /*    40 */   187,   74,   74,   74,   74,   74,   74,   74,  192,  194,
 /*    50 */   116,  100,  137,  102,  151,  153,  134,  157,  175,  177,
 /*    60 */   165,  189,  167,  205,  195,  216,  207,  218,  210,  213,
 /*    70 */   213,  213,
};
#define YY_REDUCE_USE_DFLT (-25)
#define YY_REDUCE_MAX 30
static const short yy_reduce_ofst[] = {
 /*     0 */    39,   70,   97,  108,  114,  122,  128,  135,  142,  148,
 /*    10 */   156,  162,   87,  115,  129,  136,  149,  163,  169,  170,
 /*    20 */   176,   21,  -12,  -24,  -21,   56,   63,   -2,   62,   65,
 /*    30 */    71,
};
static const YYACTIONTYPE yy_default[] = {
 /*     0 */   153,  153,  153,  153,  153,  153,  153,  153,  153,  153,
 /*    10 */   153,  153,  153,  153,  153,  153,  153,  153,  153,  153,
 /*    20 */   153,  144,  105,  153,  153,  153,  153,  153,  153,  153,
 /*    30 */   153,  153,  153,  153,  153,  118,  148,  147,  153,  153,
 /*    40 */   153,  111,  112,  113,  114,  115,  116,  117,  153,  153,
 /*    50 */   125,  153,  129,  153,  119,  123,  107,  153,  153,  153,
 /*    60 */   153,  153,  153,  153,  153,  153,  153,  153,  153,  149,
 /*    70 */   150,  151,  104,  126,  131,  132,  127,  128,  130,  133,
 /*    80 */   134,  120,  121,  122,  124,  106,  108,  109,  110,  135,
 /*    90 */   136,  137,  138,  143,  144,  139,  140,  141,  142,  145,
 /*   100 */   146,  152,
};
#define YY_SZ_ACTTAB (int)(sizeof(yy_action)/sizeof(yy_action[0]))

/* The next table maps tokens into fallback tokens.  If a construct
** like the following:
** 
**      %fallback ID X Y Z.
**
** appears in the grammar, then ID becomes a fallback token for X, Y,
** and Z.  Whenever one of the tokens X, Y, or Z is input to the parser
** but it does not parse, the type of the token is changed to ID and
** the parse is retried before an error is thrown.
*/
#ifdef YYFALLBACK
static const YYCODETYPE yyFallback[] = {
};
#endif /* YYFALLBACK */

/* The following structure represents a single element of the
** parser's stack.  Information stored includes:
**
**   +  The state number for the parser at this level of the stack.
**
**   +  The value of the token stored at this level of the stack.
**      (In other words, the "major" token.)
**
**   +  The semantic value stored at this level of the stack.  This is
**      the information used by the action routines in the grammar.
**      It is sometimes called the "minor" token.
*/
struct yyStackEntry {
  YYACTIONTYPE stateno;  /* The state-number */
  YYCODETYPE major;      /* The major token value.  This is the code
                         ** number for the token at this stack level */
  YYMINORTYPE minor;     /* The user-supplied minor token value.  This
                         ** is the value of the token  */
};
typedef struct yyStackEntry yyStackEntry;

/* The state of the parser is completely contained in an instance of
** the following structure */
struct yyParser {
  int yyidx;                    /* Index of top element in stack */
#ifdef YYTRACKMAXSTACKDEPTH
  int yyidxMax;                 /* Maximum value of yyidx */
#endif
  int yyerrcnt;                 /* Shifts left before out of the error */
  atomselParseARG_SDECL                /* A place to hold %extra_argument */
#if YYSTACKDEPTH<=0
  int yystksz;                  /* Current side of the stack */
  yyStackEntry *yystack;        /* The parser's stack */
#else
  yyStackEntry yystack[YYSTACKDEPTH];  /* The parser's stack */
#endif
};
typedef struct yyParser yyParser;

#ifndef NDEBUG
#include <stdio.h>
static FILE *yyTraceFILE = 0;
static char *yyTracePrompt = 0;
#endif /* NDEBUG */

#ifndef NDEBUG
/* 
** Turn parser tracing on by giving a stream to which to write the trace
** and a prompt to preface each trace message.  Tracing is turned off
** by making either argument NULL 
**
** Inputs:
** <ul>
** <li> A FILE* to which trace output should be written.
**      If NULL, then tracing is turned off.
** <li> A prefix string written at the beginning of every
**      line of trace output.  If NULL, then tracing is
**      turned off.
** </ul>
**
** Outputs:
** None.
*/
void atomselParseTrace(FILE *TraceFILE, char *zTracePrompt){
  yyTraceFILE = TraceFILE;
  yyTracePrompt = zTracePrompt;
  if( yyTraceFILE==0 ) yyTracePrompt = 0;
  else if( yyTracePrompt==0 ) yyTraceFILE = 0;
}
#endif /* NDEBUG */

#ifndef NDEBUG
/* For tracing shifts, the names of all terminals and nonterminals
** are required.  The following table supplies these names */
static const char *const yyTokenName[] = { 
  "$",             "OF",            "AS",            "TO",          
  "OR",            "AND",           "ADD",           "SUB",         
  "MUL",           "DIV",           "MOD",           "EXP",         
  "NEG",           "NOT",           "VAL",           "KEY",         
  "LPAREN",        "RPAREN",        "MACRO",         "WITHIN",      
  "EXWITHIN",      "PBWITHIN",      "NEAREST",       "INT",         
  "WITHINBONDS",   "PBNEAREST",     "SAME",          "CMP",         
  "REGEX",         "FLT",           "PLUS",          "FUNC",        
  "error",         "expr",          "pexpr",         "nexpr",       
  "binexpr",       "selection",     "num",           "integer",     
  "flt",           "list",          "intlist",       "fltlist",     
  "strlist",       "input",       
};
#endif /* NDEBUG */

#ifndef NDEBUG
/* For tracing reduce actions, the names of all rules are required.
*/
static const char *const yyRuleName[] = {
 /*   0 */ "input ::= selection",
 /*   1 */ "input ::=",
 /*   2 */ "selection ::= VAL",
 /*   3 */ "selection ::= KEY list",
 /*   4 */ "selection ::= selection AND selection",
 /*   5 */ "selection ::= selection OR selection",
 /*   6 */ "selection ::= LPAREN selection RPAREN",
 /*   7 */ "selection ::= MACRO",
 /*   8 */ "selection ::= NOT selection",
 /*   9 */ "selection ::= WITHIN num OF selection",
 /*  10 */ "selection ::= EXWITHIN num OF selection",
 /*  11 */ "selection ::= PBWITHIN num OF selection",
 /*  12 */ "selection ::= NEAREST INT TO selection",
 /*  13 */ "selection ::= WITHINBONDS INT OF selection",
 /*  14 */ "selection ::= PBNEAREST INT TO selection",
 /*  15 */ "selection ::= SAME KEY AS selection",
 /*  16 */ "selection ::= expr CMP expr",
 /*  17 */ "list ::= integer",
 /*  18 */ "list ::= integer TO integer",
 /*  19 */ "list ::= VAL",
 /*  20 */ "list ::= REGEX",
 /*  21 */ "list ::= flt",
 /*  22 */ "list ::= flt TO flt",
 /*  23 */ "list ::= list integer",
 /*  24 */ "list ::= list integer TO integer",
 /*  25 */ "list ::= list VAL",
 /*  26 */ "list ::= list REGEX",
 /*  27 */ "list ::= list flt",
 /*  28 */ "list ::= list flt TO flt",
 /*  29 */ "integer ::= INT",
 /*  30 */ "integer ::= NEG INT",
 /*  31 */ "flt ::= FLT",
 /*  32 */ "flt ::= NEG FLT",
 /*  33 */ "num ::= INT",
 /*  34 */ "num ::= FLT",
 /*  35 */ "num ::= NEG INT",
 /*  36 */ "num ::= NEG FLT",
 /*  37 */ "num ::= PLUS INT",
 /*  38 */ "num ::= PLUS FLT",
 /*  39 */ "expr ::= num",
 /*  40 */ "expr ::= nexpr",
 /*  41 */ "expr ::= NEG nexpr",
 /*  42 */ "nexpr ::= KEY",
 /*  43 */ "nexpr ::= LPAREN expr RPAREN",
 /*  44 */ "nexpr ::= FUNC LPAREN expr RPAREN",
 /*  45 */ "expr ::= expr ADD expr",
 /*  46 */ "expr ::= expr SUB expr",
 /*  47 */ "expr ::= expr MUL expr",
 /*  48 */ "expr ::= expr DIV expr",
 /*  49 */ "expr ::= expr MOD expr",
 /*  50 */ "expr ::= expr EXP expr",
};
#endif /* NDEBUG */


#if YYSTACKDEPTH<=0
/*
** Try to increase the size of the parser stack.
*/
static void yyGrowStack(yyParser *p){
  int newSize;
  yyStackEntry *pNew;

  newSize = p->yystksz*2 + 100;
  pNew = realloc(p->yystack, newSize*sizeof(pNew[0]));
  if( pNew ){
    p->yystack = pNew;
    p->yystksz = newSize;
#ifndef NDEBUG
    if( yyTraceFILE ){
      fprintf(yyTraceFILE,"%sStack grows to %d entries!\n",
              yyTracePrompt, p->yystksz);
    }
#endif
  }
}
#endif

/* 
** This function allocates a new parser.
** The only argument is a pointer to a function which works like
** malloc.
**
** Inputs:
** A pointer to the function used to allocate memory.
**
** Outputs:
** A pointer to a parser.  This pointer is used in subsequent calls
** to atomselParse and atomselParseFree.
*/
void *atomselParseAlloc(void *(*mallocProc)(size_t)){
  yyParser *pParser;
  pParser = (yyParser*)(*mallocProc)( (size_t)sizeof(yyParser) );
  if( pParser ){
    pParser->yyidx = -1;
#ifdef YYTRACKMAXSTACKDEPTH
    pParser->yyidxMax = 0;
#endif
#if YYSTACKDEPTH<=0
    pParser->yystack = NULL;
    pParser->yystksz = 0;
    yyGrowStack(pParser);
#endif
  }
  return pParser;
}

/* The following function deletes the value associated with a
** symbol.  The symbol can be either a terminal or nonterminal.
** "yymajor" is the symbol code, and "yypminor" is a pointer to
** the value.
*/
static void yy_destructor(
  yyParser *yypParser,    /* The parser */
  YYCODETYPE yymajor,     /* Type code for object to destroy */
  YYMINORTYPE *yypminor   /* The object to be destroyed */
){
  atomselParseARG_FETCH;
  switch( yymajor ){
    /* Here is inserted the actions which take place when a
    ** terminal or non-terminal is destroyed.  This can happen
    ** when the symbol is popped from the stack during a
    ** reduce or during error processing or when a parser is 
    ** being destroyed before it is finished parsing.
    **
    ** Note: during a reduce, the only symbols destroyed are those
    ** which appear on the RHS of the rule, but which are not used
    ** inside the C code.
    */
    case 33: /* expr */
    case 34: /* pexpr */
    case 35: /* nexpr */
{
#line 38 "atomsel.y"
 (void)query; delete (yypminor->yy19); 
#line 484 "atomsel.c"
}
      break;
    case 37: /* selection */
{
#line 41 "atomsel.y"
 delete (yypminor->yy16); 
#line 491 "atomsel.c"
}
      break;
    case 41: /* list */
    case 42: /* intlist */
    case 43: /* fltlist */
    case 44: /* strlist */
{
#line 42 "atomsel.y"
 delete (yypminor->yy14); 
#line 501 "atomsel.c"
}
      break;
    default:  break;   /* If no destructor action specified: do nothing */
  }
}

/*
** Pop the parser's stack once.
**
** If there is a destructor routine associated with the token which
** is popped from the stack, then call it.
**
** Return the major token number for the symbol popped.
*/
static int yy_pop_parser_stack(yyParser *pParser){
  YYCODETYPE yymajor;
  yyStackEntry *yytos = &pParser->yystack[pParser->yyidx];

  if( pParser->yyidx<0 ) return 0;
#ifndef NDEBUG
  if( yyTraceFILE && pParser->yyidx>=0 ){
    fprintf(yyTraceFILE,"%sPopping %s\n",
      yyTracePrompt,
      yyTokenName[yytos->major]);
  }
#endif
  yymajor = yytos->major;
  yy_destructor(pParser, yymajor, &yytos->minor);
  pParser->yyidx--;
  return yymajor;
}

/* 
** Deallocate and destroy a parser.  Destructors are all called for
** all stack elements before shutting the parser down.
**
** Inputs:
** <ul>
** <li>  A pointer to the parser.  This should be a pointer
**       obtained from atomselParseAlloc.
** <li>  A pointer to a function used to reclaim memory obtained
**       from malloc.
** </ul>
*/
void atomselParseFree(
  void *p,                    /* The parser to be deleted */
  void (*freeProc)(void*)     /* Function used to reclaim memory */
){
  yyParser *pParser = (yyParser*)p;
  if( pParser==0 ) return;
  while( pParser->yyidx>=0 ) yy_pop_parser_stack(pParser);
#if YYSTACKDEPTH<=0
  free(pParser->yystack);
#endif
  (*freeProc)((void*)pParser);
}

/*
** Return the peak depth of the stack for a parser.
*/
#ifdef YYTRACKMAXSTACKDEPTH
int atomselParseStackPeak(void *p){
  yyParser *pParser = (yyParser*)p;
  return pParser->yyidxMax;
}
#endif

/*
** Find the appropriate action for a parser given the terminal
** look-ahead token iLookAhead.
**
** If the look-ahead token is YYNOCODE, then check to see if the action is
** independent of the look-ahead.  If it is, return the action, otherwise
** return YY_NO_ACTION.
*/
static int yy_find_shift_action(
  yyParser *pParser,        /* The parser */
  YYCODETYPE iLookAhead     /* The look-ahead token */
){
  int i;
  int stateno = pParser->yystack[pParser->yyidx].stateno;
 
  if( stateno>YY_SHIFT_MAX || (i = yy_shift_ofst[stateno])==YY_SHIFT_USE_DFLT ){
    return yy_default[stateno];
  }
  assert( iLookAhead!=YYNOCODE );
  i += iLookAhead;
  if( i<0 || i>=YY_SZ_ACTTAB || yy_lookahead[i]!=iLookAhead ){
    if( iLookAhead>0 ){
#ifdef YYFALLBACK
      YYCODETYPE iFallback;            /* Fallback token */
      if( iLookAhead<sizeof(yyFallback)/sizeof(yyFallback[0])
             && (iFallback = yyFallback[iLookAhead])!=0 ){
#ifndef NDEBUG
        if( yyTraceFILE ){
          fprintf(yyTraceFILE, "%sFALLBACK %s => %s\n",
             yyTracePrompt, yyTokenName[iLookAhead], yyTokenName[iFallback]);
        }
#endif
        return yy_find_shift_action(pParser, iFallback);
      }
#endif
#ifdef YYWILDCARD
      {
        int j = i - iLookAhead + YYWILDCARD;
        if( j>=0 && j<YY_SZ_ACTTAB && yy_lookahead[j]==YYWILDCARD ){
#ifndef NDEBUG
          if( yyTraceFILE ){
            fprintf(yyTraceFILE, "%sWILDCARD %s => %s\n",
               yyTracePrompt, yyTokenName[iLookAhead], yyTokenName[YYWILDCARD]);
          }
#endif /* NDEBUG */
          return yy_action[j];
        }
      }
#endif /* YYWILDCARD */
    }
    return yy_default[stateno];
  }else{
    return yy_action[i];
  }
}

/*
** Find the appropriate action for a parser given the non-terminal
** look-ahead token iLookAhead.
**
** If the look-ahead token is YYNOCODE, then check to see if the action is
** independent of the look-ahead.  If it is, return the action, otherwise
** return YY_NO_ACTION.
*/
static int yy_find_reduce_action(
  int stateno,              /* Current state number */
  YYCODETYPE iLookAhead     /* The look-ahead token */
){
  int i;
#ifdef YYERRORSYMBOL
  if( stateno>YY_REDUCE_MAX ){
    return yy_default[stateno];
  }
#else
  assert( stateno<=YY_REDUCE_MAX );
#endif
  i = yy_reduce_ofst[stateno];
  assert( i!=YY_REDUCE_USE_DFLT );
  assert( iLookAhead!=YYNOCODE );
  i += iLookAhead;
#ifdef YYERRORSYMBOL
  if( i<0 || i>=YY_SZ_ACTTAB || yy_lookahead[i]!=iLookAhead ){
    return yy_default[stateno];
  }
#else
  assert( i>=0 && i<YY_SZ_ACTTAB );
  assert( yy_lookahead[i]==iLookAhead );
#endif
  return yy_action[i];
}

/*
** The following routine is called if the stack overflows.
*/
static void yyStackOverflow(yyParser *yypParser, YYMINORTYPE *yypMinor){
   atomselParseARG_FETCH;
   yypParser->yyidx--;
#ifndef NDEBUG
   if( yyTraceFILE ){
     fprintf(yyTraceFILE,"%sStack Overflow!\n",yyTracePrompt);
   }
#endif
   while( yypParser->yyidx>=0 ) yy_pop_parser_stack(yypParser);
   /* Here code is inserted which will execute if the parser
   ** stack every overflows */
   atomselParseARG_STORE; /* Suppress warning about unused %extra_argument var */
}

/*
** Perform a shift action.
*/
static void yy_shift(
  yyParser *yypParser,          /* The parser to be shifted */
  int yyNewState,               /* The new state to shift in */
  int yyMajor,                  /* The major token to shift in */
  YYMINORTYPE *yypMinor         /* Pointer to the minor token to shift in */
){
  yyStackEntry *yytos;
  yypParser->yyidx++;
#ifdef YYTRACKMAXSTACKDEPTH
  if( yypParser->yyidx>yypParser->yyidxMax ){
    yypParser->yyidxMax = yypParser->yyidx;
  }
#endif
#if YYSTACKDEPTH>0 
  if( yypParser->yyidx>=YYSTACKDEPTH ){
    yyStackOverflow(yypParser, yypMinor);
    return;
  }
#else
  if( yypParser->yyidx>=yypParser->yystksz ){
    yyGrowStack(yypParser);
    if( yypParser->yyidx>=yypParser->yystksz ){
      yyStackOverflow(yypParser, yypMinor);
      return;
    }
  }
#endif
  yytos = &yypParser->yystack[yypParser->yyidx];
  yytos->stateno = (YYACTIONTYPE)yyNewState;
  yytos->major = (YYCODETYPE)yyMajor;
  yytos->minor = *yypMinor;
#ifndef NDEBUG
  if( yyTraceFILE && yypParser->yyidx>0 ){
    int i;
    fprintf(yyTraceFILE,"%sShift %d\n",yyTracePrompt,yyNewState);
    fprintf(yyTraceFILE,"%sStack:",yyTracePrompt);
    for(i=1; i<=yypParser->yyidx; i++)
      fprintf(yyTraceFILE," %s",yyTokenName[yypParser->yystack[i].major]);
    fprintf(yyTraceFILE,"\n");
  }
#endif
}

/* The following table contains information about every rule that
** is used during the reduce.
*/
static const struct {
  YYCODETYPE lhs;         /* Symbol on the left-hand side of the rule */
  unsigned char nrhs;     /* Number of right-hand side symbols in the rule */
} yyRuleInfo[] = {
  { 45, 1 },
  { 45, 0 },
  { 37, 1 },
  { 37, 2 },
  { 37, 3 },
  { 37, 3 },
  { 37, 3 },
  { 37, 1 },
  { 37, 2 },
  { 37, 4 },
  { 37, 4 },
  { 37, 4 },
  { 37, 4 },
  { 37, 4 },
  { 37, 4 },
  { 37, 4 },
  { 37, 3 },
  { 41, 1 },
  { 41, 3 },
  { 41, 1 },
  { 41, 1 },
  { 41, 1 },
  { 41, 3 },
  { 41, 2 },
  { 41, 4 },
  { 41, 2 },
  { 41, 2 },
  { 41, 2 },
  { 41, 4 },
  { 39, 1 },
  { 39, 2 },
  { 40, 1 },
  { 40, 2 },
  { 38, 1 },
  { 38, 1 },
  { 38, 2 },
  { 38, 2 },
  { 38, 2 },
  { 38, 2 },
  { 33, 1 },
  { 33, 1 },
  { 33, 2 },
  { 35, 1 },
  { 35, 3 },
  { 35, 4 },
  { 33, 3 },
  { 33, 3 },
  { 33, 3 },
  { 33, 3 },
  { 33, 3 },
  { 33, 3 },
};

static void yy_accept(yyParser*);  /* Forward Declaration */

/*
** Perform a reduce action and the shift that must immediately
** follow the reduce.
*/
static void yy_reduce(
  yyParser *yypParser,         /* The parser */
  int yyruleno                 /* Number of the rule by which to reduce */
){
  int yygoto;                     /* The next state */
  int yyact;                      /* The next action */
  YYMINORTYPE yygotominor;        /* The LHS of the rule reduced */
  yyStackEntry *yymsp;            /* The top of the parser's stack */
  int yysize;                     /* Amount to pop the stack */
  atomselParseARG_FETCH;
  yymsp = &yypParser->yystack[yypParser->yyidx];
#ifndef NDEBUG
  if( yyTraceFILE && yyruleno>=0 
        && yyruleno<(int)(sizeof(yyRuleName)/sizeof(yyRuleName[0])) ){
    fprintf(yyTraceFILE, "%sReduce [%s].\n", yyTracePrompt,
      yyRuleName[yyruleno]);
  }
#endif /* NDEBUG */

  /* Silence complaints from purify about yygotominor being uninitialized
  ** in some cases when it is copied into the stack after the following
  ** switch.  yygotominor is uninitialized when a rule reduces that does
  ** not set the value of its left-hand side nonterminal.  Leaving the
  ** value of the nonterminal uninitialized is utterly harmless as long
  ** as the value is never used.  So really the only thing this code
  ** accomplishes is to quieten purify.  
  **
  ** 2007-01-16:  The wireshark project (www.wireshark.org) reports that
  ** without this code, their parser segfaults.  I'm not sure what there
  ** parser is doing to make this happen.  This is the second bug report
  ** from wireshark this week.  Clearly they are stressing Lemon in ways
  ** that it has not been previously stressed...  (SQLite ticket #2172)
  */
  /*memset(&yygotominor, 0, sizeof(yygotominor));*/
  yygotominor = yyzerominor;


  switch( yyruleno ){
  /* Beginning here are the reduction cases.  A typical example
  ** follows:
  **   case 0:
  **  #line <lineno> <grammarfile>
  **     { ... }           // User supplied code
  **  #line <lineno> <thisfile>
  **     break;
  */
      case 0: /* input ::= selection */
#line 49 "atomsel.y"
{ query->pred.reset(yymsp[0].minor.yy16); }
#line 838 "atomsel.c"
        break;
      case 2: /* selection ::= VAL */
#line 54 "atomsel.y"
{ yygotominor.yy16=new BoolPredicate(query->mol,yymsp[0].minor.yy0.str());  }
#line 843 "atomsel.c"
        break;
      case 3: /* selection ::= KEY list */
#line 55 "atomsel.y"
{ yygotominor.yy16=new KeyPredicate(query,yymsp[-1].minor.yy0.str(),yymsp[0].minor.yy14); }
#line 848 "atomsel.c"
        break;
      case 4: /* selection ::= selection AND selection */
#line 56 "atomsel.y"
{yygotominor.yy16=new AndPredicate(yymsp[-2].minor.yy16,yymsp[0].minor.yy16); }
#line 853 "atomsel.c"
        break;
      case 5: /* selection ::= selection OR selection */
#line 57 "atomsel.y"
{yygotominor.yy16=new OrPredicate(yymsp[-2].minor.yy16,yymsp[0].minor.yy16); }
#line 858 "atomsel.c"
        break;
      case 6: /* selection ::= LPAREN selection RPAREN */
#line 58 "atomsel.y"
{yygotominor.yy16=yymsp[-1].minor.yy16; }
#line 863 "atomsel.c"
        break;
      case 7: /* selection ::= MACRO */
#line 59 "atomsel.y"
{yygotominor.yy16=query->pred.release(); }
#line 868 "atomsel.c"
        break;
      case 8: /* selection ::= NOT selection */
#line 60 "atomsel.y"
{yygotominor.yy16=new NotPredicate(yymsp[0].minor.yy16); }
#line 873 "atomsel.c"
        break;
      case 9: /* selection ::= WITHIN num OF selection */
#line 61 "atomsel.y"
{yygotominor.yy16=new WithinPredicate(query->mol, query->pos, NULL, yymsp[-2].minor.yy76, false, false, yymsp[0].minor.yy16); }
#line 878 "atomsel.c"
        break;
      case 10: /* selection ::= EXWITHIN num OF selection */
#line 62 "atomsel.y"
{yygotominor.yy16=new WithinPredicate(query->mol, query->pos, NULL, yymsp[-2].minor.yy76,  true, false, yymsp[0].minor.yy16); }
#line 883 "atomsel.c"
        break;
      case 11: /* selection ::= PBWITHIN num OF selection */
#line 63 "atomsel.y"
{yygotominor.yy16=new WithinPredicate(query->mol, query->pos, query->cell, yymsp[-2].minor.yy76, false,  true, yymsp[0].minor.yy16); }
#line 888 "atomsel.c"
        break;
      case 12: /* selection ::= NEAREST INT TO selection */
#line 64 "atomsel.y"
{yygotominor.yy16=new KNearestPredicate(query->mol,query->pos, NULL, yymsp[-2].minor.yy0.ival, false, yymsp[0].minor.yy16); }
#line 893 "atomsel.c"
        break;
      case 13: /* selection ::= WITHINBONDS INT OF selection */
#line 65 "atomsel.y"
{yygotominor.yy16=new WithinBondsPredicate(query->mol,yymsp[-2].minor.yy0.ival, yymsp[0].minor.yy16); }
#line 898 "atomsel.c"
        break;
      case 14: /* selection ::= PBNEAREST INT TO selection */
#line 66 "atomsel.y"
{yygotominor.yy16=new KNearestPredicate(query->mol,query->pos,query->cell, yymsp[-2].minor.yy0.ival,  true, yymsp[0].minor.yy16); }
#line 903 "atomsel.c"
        break;
      case 15: /* selection ::= SAME KEY AS selection */
#line 67 "atomsel.y"
{yygotominor.yy16=new SamePredicate(query,yymsp[-2].minor.yy0.str(),yymsp[0].minor.yy16); }
#line 908 "atomsel.c"
        break;
      case 16: /* selection ::= expr CMP expr */
#line 68 "atomsel.y"
{yygotominor.yy16=new CmpPredicate(yymsp[-1].minor.yy0.ival,yymsp[-2].minor.yy19,yymsp[0].minor.yy19);}
#line 913 "atomsel.c"
        break;
      case 17: /* list ::= integer */
#line 70 "atomsel.y"
{ yygotominor.yy14=new Valist; yygotominor.yy14->add(yymsp[0].minor.yy52); }
#line 918 "atomsel.c"
        break;
      case 18: /* list ::= integer TO integer */
#line 71 "atomsel.y"
{ yygotominor.yy14=new Valist; yygotominor.yy14->add(yymsp[-2].minor.yy52,yymsp[0].minor.yy52); }
#line 923 "atomsel.c"
        break;
      case 19: /* list ::= VAL */
#line 72 "atomsel.y"
{ yygotominor.yy14=new Valist; yygotominor.yy14->add(yymsp[0].minor.yy0.str()); }
#line 928 "atomsel.c"
        break;
      case 20: /* list ::= REGEX */
#line 73 "atomsel.y"
{ yygotominor.yy14=new Valist; yygotominor.yy14->add_regex(yymsp[0].minor.yy0.str()); }
#line 933 "atomsel.c"
        break;
      case 21: /* list ::= flt */
#line 74 "atomsel.y"
{ yygotominor.yy14=new Valist; yygotominor.yy14->add(yymsp[0].minor.yy76); }
#line 938 "atomsel.c"
        break;
      case 22: /* list ::= flt TO flt */
#line 75 "atomsel.y"
{ yygotominor.yy14=new Valist; yygotominor.yy14->add(yymsp[-2].minor.yy76,yymsp[0].minor.yy76); }
#line 943 "atomsel.c"
        break;
      case 23: /* list ::= list integer */
#line 77 "atomsel.y"
{ yygotominor.yy14=yymsp[-1].minor.yy14; yygotominor.yy14->add(yymsp[0].minor.yy52); }
#line 948 "atomsel.c"
        break;
      case 24: /* list ::= list integer TO integer */
#line 78 "atomsel.y"
{ yygotominor.yy14=yymsp[-3].minor.yy14; yygotominor.yy14->add(yymsp[-2].minor.yy52,yymsp[0].minor.yy52); }
#line 953 "atomsel.c"
        break;
      case 25: /* list ::= list VAL */
#line 79 "atomsel.y"
{ yygotominor.yy14=yymsp[-1].minor.yy14; yygotominor.yy14->add(yymsp[0].minor.yy0.str()); }
#line 958 "atomsel.c"
        break;
      case 26: /* list ::= list REGEX */
#line 80 "atomsel.y"
{ yygotominor.yy14=yymsp[-1].minor.yy14; yygotominor.yy14->add_regex(yymsp[0].minor.yy0.str()); }
#line 963 "atomsel.c"
        break;
      case 27: /* list ::= list flt */
#line 81 "atomsel.y"
{ yygotominor.yy14=yymsp[-1].minor.yy14; yygotominor.yy14->add(yymsp[0].minor.yy76); }
#line 968 "atomsel.c"
        break;
      case 28: /* list ::= list flt TO flt */
#line 82 "atomsel.y"
{ yygotominor.yy14=yymsp[-3].minor.yy14; yygotominor.yy14->add(yymsp[-2].minor.yy76,yymsp[0].minor.yy76); }
#line 973 "atomsel.c"
        break;
      case 29: /* integer ::= INT */
#line 84 "atomsel.y"
{ yygotominor.yy52= yymsp[0].minor.yy0.ival; }
#line 978 "atomsel.c"
        break;
      case 30: /* integer ::= NEG INT */
#line 85 "atomsel.y"
{ yygotominor.yy52=-yymsp[0].minor.yy0.ival; }
#line 983 "atomsel.c"
        break;
      case 31: /* flt ::= FLT */
#line 87 "atomsel.y"
{ yygotominor.yy76= yymsp[0].minor.yy0.fval; }
#line 988 "atomsel.c"
        break;
      case 32: /* flt ::= NEG FLT */
      case 36: /* num ::= NEG FLT */ yytestcase(yyruleno==36);
#line 88 "atomsel.y"
{ yygotominor.yy76=-yymsp[0].minor.yy0.fval; }
#line 994 "atomsel.c"
        break;
      case 33: /* num ::= INT */
#line 90 "atomsel.y"
{ yygotominor.yy76=yymsp[0].minor.yy0.ival;  }
#line 999 "atomsel.c"
        break;
      case 34: /* num ::= FLT */
#line 91 "atomsel.y"
{ yygotominor.yy76=yymsp[0].minor.yy0.fval;  }
#line 1004 "atomsel.c"
        break;
      case 35: /* num ::= NEG INT */
#line 92 "atomsel.y"
{ yygotominor.yy76=-yymsp[0].minor.yy0.ival; }
#line 1009 "atomsel.c"
        break;
      case 37: /* num ::= PLUS INT */
#line 94 "atomsel.y"
{ yygotominor.yy76=yymsp[0].minor.yy0.ival; }
#line 1014 "atomsel.c"
        break;
      case 38: /* num ::= PLUS FLT */
#line 95 "atomsel.y"
{ yygotominor.yy76=yymsp[0].minor.yy0.fval; }
#line 1019 "atomsel.c"
        break;
      case 39: /* expr ::= num */
#line 97 "atomsel.y"
{ yygotominor.yy19=new LitExpr(yymsp[0].minor.yy76); }
#line 1024 "atomsel.c"
        break;
      case 40: /* expr ::= nexpr */
#line 98 "atomsel.y"
{ yygotominor.yy19=yymsp[0].minor.yy19; }
#line 1029 "atomsel.c"
        break;
      case 41: /* expr ::= NEG nexpr */
#line 99 "atomsel.y"
{ yygotominor.yy19=new NegExpr(yymsp[0].minor.yy19); }
#line 1034 "atomsel.c"
        break;
      case 42: /* nexpr ::= KEY */
#line 101 "atomsel.y"
{ yygotominor.yy19=new KeyExpr(query,yymsp[0].minor.yy0.str()); }
#line 1039 "atomsel.c"
        break;
      case 43: /* nexpr ::= LPAREN expr RPAREN */
#line 102 "atomsel.y"
{ yygotominor.yy19=yymsp[-1].minor.yy19; }
#line 1044 "atomsel.c"
        break;
      case 44: /* nexpr ::= FUNC LPAREN expr RPAREN */
#line 103 "atomsel.y"
{ yygotominor.yy19=new FuncExpr(yymsp[-3].minor.yy0.func, yymsp[-1].minor.yy19); }
#line 1049 "atomsel.c"
        break;
      case 45: /* expr ::= expr ADD expr */
#line 105 "atomsel.y"
{yygotominor.yy19=new BinExpr(ADD,yymsp[-2].minor.yy19,yymsp[0].minor.yy19); }
#line 1054 "atomsel.c"
        break;
      case 46: /* expr ::= expr SUB expr */
#line 106 "atomsel.y"
{yygotominor.yy19=new BinExpr(SUB,yymsp[-2].minor.yy19,yymsp[0].minor.yy19); }
#line 1059 "atomsel.c"
        break;
      case 47: /* expr ::= expr MUL expr */
#line 107 "atomsel.y"
{yygotominor.yy19=new BinExpr(MUL,yymsp[-2].minor.yy19,yymsp[0].minor.yy19); }
#line 1064 "atomsel.c"
        break;
      case 48: /* expr ::= expr DIV expr */
#line 108 "atomsel.y"
{yygotominor.yy19=new BinExpr(DIV,yymsp[-2].minor.yy19,yymsp[0].minor.yy19); }
#line 1069 "atomsel.c"
        break;
      case 49: /* expr ::= expr MOD expr */
#line 109 "atomsel.y"
{yygotominor.yy19=new BinExpr(MOD,yymsp[-2].minor.yy19,yymsp[0].minor.yy19); }
#line 1074 "atomsel.c"
        break;
      case 50: /* expr ::= expr EXP expr */
#line 110 "atomsel.y"
{yygotominor.yy19=new BinExpr(EXP,yymsp[-2].minor.yy19,yymsp[0].minor.yy19); }
#line 1079 "atomsel.c"
        break;
      default:
      /* (1) input ::= */ yytestcase(yyruleno==1);
        break;
  };
  yygoto = yyRuleInfo[yyruleno].lhs;
  yysize = yyRuleInfo[yyruleno].nrhs;
  yypParser->yyidx -= yysize;
  yyact = yy_find_reduce_action(yymsp[-yysize].stateno,(YYCODETYPE)yygoto);
  if( yyact < YYNSTATE ){
#ifdef NDEBUG
    /* If we are not debugging and the reduce action popped at least
    ** one element off the stack, then we can push the new element back
    ** onto the stack here, and skip the stack overflow test in yy_shift().
    ** That gives a significant speed improvement. */
    if( yysize ){
      yypParser->yyidx++;
      yymsp -= yysize-1;
      yymsp->stateno = (YYACTIONTYPE)yyact;
      yymsp->major = (YYCODETYPE)yygoto;
      yymsp->minor = yygotominor;
    }else
#endif
    {
      yy_shift(yypParser,yyact,yygoto,&yygotominor);
    }
  }else{
    assert( yyact == YYNSTATE + YYNRULE + 1 );
    yy_accept(yypParser);
  }
}

/*
** The following code executes when the parse fails
*/
#ifndef YYNOERRORRECOVERY
static void yy_parse_failed(
  yyParser *yypParser           /* The parser */
){
  atomselParseARG_FETCH;
#ifndef NDEBUG
  if( yyTraceFILE ){
    fprintf(yyTraceFILE,"%sFail!\n",yyTracePrompt);
  }
#endif
  while( yypParser->yyidx>=0 ) yy_pop_parser_stack(yypParser);
  /* Here code is inserted which will be executed whenever the
  ** parser fails */
  atomselParseARG_STORE; /* Suppress warning about unused %extra_argument variable */
}
#endif /* YYNOERRORRECOVERY */

/*
** The following code executes when a syntax error first occurs.
*/
static void yy_syntax_error(
  yyParser *yypParser,           /* The parser */
  int yymajor,                   /* The major type of the error token */
  YYMINORTYPE yyminor            /* The minor type of the error token */
){
  atomselParseARG_FETCH;
#define TOKEN (yyminor.yy0)
#line 47 "atomsel.y"
 query->mol = nullptr; 
#line 1144 "atomsel.c"
  atomselParseARG_STORE; /* Suppress warning about unused %extra_argument variable */
}

/*
** The following is executed when the parser accepts
*/
static void yy_accept(
  yyParser *yypParser           /* The parser */
){
  atomselParseARG_FETCH;
#ifndef NDEBUG
  if( yyTraceFILE ){
    fprintf(yyTraceFILE,"%sAccept!\n",yyTracePrompt);
  }
#endif
  while( yypParser->yyidx>=0 ) yy_pop_parser_stack(yypParser);
  /* Here code is inserted which will be executed whenever the
  ** parser accepts */
  atomselParseARG_STORE; /* Suppress warning about unused %extra_argument variable */
}

/* The main parser program.
** The first argument is a pointer to a structure obtained from
** "atomselParseAlloc" which describes the current state of the parser.
** The second argument is the major token number.  The third is
** the minor token.  The fourth optional argument is whatever the
** user wants (and specified in the grammar) and is available for
** use by the action routines.
**
** Inputs:
** <ul>
** <li> A pointer to the parser (an opaque structure.)
** <li> The major token number.
** <li> The minor token number.
** <li> An option argument of a grammar-specified type.
** </ul>
**
** Outputs:
** None.
*/
void atomselParse(
  void *yyp,                   /* The parser */
  int yymajor,                 /* The major token code number */
  atomselParseTOKENTYPE yyminor       /* The value for the token */
  atomselParseARG_PDECL               /* Optional %extra_argument parameter */
){
  YYMINORTYPE yyminorunion;
  int yyact;            /* The parser action. */
  int yyendofinput;     /* True if we are at the end of input */
#ifdef YYERRORSYMBOL
  int yyerrorhit = 0;   /* True if yymajor has invoked an error */
#endif
  yyParser *yypParser;  /* The parser */

  /* (re)initialize the parser, if necessary */
  yypParser = (yyParser*)yyp;
  if( yypParser->yyidx<0 ){
#if YYSTACKDEPTH<=0
    if( yypParser->yystksz <=0 ){
      /*memset(&yyminorunion, 0, sizeof(yyminorunion));*/
      yyminorunion = yyzerominor;
      yyStackOverflow(yypParser, &yyminorunion);
      return;
    }
#endif
    yypParser->yyidx = 0;
    yypParser->yyerrcnt = -1;
    yypParser->yystack[0].stateno = 0;
    yypParser->yystack[0].major = 0;
  }
  yyminorunion.yy0 = yyminor;
  yyendofinput = (yymajor==0) ? true : false;
  atomselParseARG_STORE;

#ifndef NDEBUG
  if( yyTraceFILE ){
    fprintf(yyTraceFILE,"%sInput %s\n",yyTracePrompt,yyTokenName[yymajor]);
  }
#endif

  do{
    yyact = yy_find_shift_action(yypParser,(YYCODETYPE)yymajor);
    if( yyact<YYNSTATE ){
      assert( !yyendofinput );  /* Impossible to shift the $ token */
      yy_shift(yypParser,yyact,yymajor,&yyminorunion);
      yypParser->yyerrcnt--;
      yymajor = YYNOCODE;
    }else if( yyact < YYNSTATE + YYNRULE ){
      yy_reduce(yypParser,yyact-YYNSTATE);
    }else{
      assert( yyact == YY_ERROR_ACTION );
#ifdef YYERRORSYMBOL
      int yymx;
#endif
#ifndef NDEBUG
      if( yyTraceFILE ){
        fprintf(yyTraceFILE,"%sSyntax Error!\n",yyTracePrompt);
      }
#endif
#ifdef YYERRORSYMBOL
      /* A syntax error has occurred.
      ** The response to an error depends upon whether or not the
      ** grammar defines an error token "ERROR".  
      **
      ** This is what we do if the grammar does define ERROR:
      **
      **  * Call the %syntax_error function.
      **
      **  * Begin popping the stack until we enter a state where
      **    it is legal to shift the error symbol, then shift
      **    the error symbol.
      **
      **  * Set the error count to three.
      **
      **  * Begin accepting and shifting new tokens.  No new error
      **    processing will occur until three tokens have been
      **    shifted successfully.
      **
      */
      if( yypParser->yyerrcnt<0 ){
        yy_syntax_error(yypParser,yymajor,yyminorunion);
      }
      yymx = yypParser->yystack[yypParser->yyidx].major;
      if( yymx==YYERRORSYMBOL || yyerrorhit ){
#ifndef NDEBUG
        if( yyTraceFILE ){
          fprintf(yyTraceFILE,"%sDiscard input token %s\n",
             yyTracePrompt,yyTokenName[yymajor]);
        }
#endif
        yy_destructor(yypParser, (YYCODETYPE)yymajor,&yyminorunion);
        yymajor = YYNOCODE;
      }else{
         while(
          yypParser->yyidx >= 0 &&
          yymx != YYERRORSYMBOL &&
          (yyact = yy_find_reduce_action(
                        yypParser->yystack[yypParser->yyidx].stateno,
                        YYERRORSYMBOL)) >= YYNSTATE
        ){
          yy_pop_parser_stack(yypParser);
        }
        if( yypParser->yyidx < 0 || yymajor==0 ){
          yy_destructor(yypParser,(YYCODETYPE)yymajor,&yyminorunion);
          yy_parse_failed(yypParser);
          yymajor = YYNOCODE;
        }else if( yymx!=YYERRORSYMBOL ){
          YYMINORTYPE u2;
          u2.YYERRSYMDT = 0;
          yy_shift(yypParser,yyact,YYERRORSYMBOL,&u2);
        }
      }
      yypParser->yyerrcnt = 3;
      yyerrorhit = 1;
#elif defined(YYNOERRORRECOVERY)
      /* If the YYNOERRORRECOVERY macro is defined, then do not attempt to
      ** do any kind of error recovery.  Instead, simply invoke the syntax
      ** error routine and continue going as if nothing had happened.
      **
      ** Applications can set this macro (for example inside %include) if
      ** they intend to abandon the parse upon the first syntax error seen.
      */
      yy_syntax_error(yypParser,yymajor,yyminorunion);
      yy_destructor(yypParser,(YYCODETYPE)yymajor,&yyminorunion);
      yymajor = YYNOCODE;
      
#else  /* YYERRORSYMBOL is not defined */
      /* This is what we do if the grammar does not define ERROR:
      **
      **  * Report an error message, and throw away the input token.
      **
      **  * If the input token is $, then fail the parse.
      **
      ** As before, subsequent error messages are suppressed until
      ** three input tokens have been successfully shifted.
      */
      if( yypParser->yyerrcnt<=0 ){
        yy_syntax_error(yypParser,yymajor,yyminorunion);
      }
      yypParser->yyerrcnt = 3;
      yy_destructor(yypParser,(YYCODETYPE)yymajor,&yyminorunion);
      if( yyendofinput ){
        yy_parse_failed(yypParser);
      }
      yymajor = YYNOCODE;
#endif
    }
  }while( yymajor!=YYNOCODE && yypParser->yyidx>=0 );
  return;
}
