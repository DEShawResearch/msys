/* A Bison parser, made by GNU Bison 2.5.  */

/* Bison implementation for Yacc-like parsers in C
   
      Copyright (C) 1984, 1989-1990, 2000-2011 Free Software Foundation, Inc.
   
   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.
   
   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
   
   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.  */

/* As a special exception, you may create a larger work that contains
   part or all of the Bison parser skeleton and distribute that work
   under terms of your choice, so long as that work isn't itself a
   parser generator using the skeleton or a modified version thereof
   as a parser skeleton.  Alternatively, if you modify or redistribute
   the parser skeleton itself, you may (at your option) remove this
   special exception, which will cause the skeleton and the resulting
   Bison output files to be licensed under the GNU General Public
   License without this special exception.
   
   This special exception was added by the Free Software Foundation in
   version 2.2 of Bison.  */

/* C LALR(1) parser skeleton written by Richard Stallman, by
   simplifying the original so-called "semantic" parser.  */

/* All symbols defined below should begin with yy or YY, to avoid
   infringing on user name space.  This should be done even for local
   variables, as they might otherwise be expanded by user macros.
   There are some unavoidable exceptions within include files to
   define necessary library symbols; they are noted "INFRINGES ON
   USER NAME SPACE" below.  */

/* Identify Bison output.  */
#define YYBISON 1

/* Bison version.  */
#define YYBISON_VERSION "2.5"

/* Skeleton name.  */
#define YYSKELETON_NAME "yacc.c"

/* Pure parsers.  */
#define YYPURE 1

/* Push parsers.  */
#define YYPUSH 0

/* Pull parsers.  */
#define YYPULL 1

/* Using locations.  */
#define YYLSP_NEEDED 1

/* Substitute the variable and function names.  */
#define yyparse         vmd_parse
#define yylex           vmd_lex
#define yyerror         vmd_error
#define yylval          vmd_lval
#define yychar          vmd_char
#define yydebug         vmd_debug
#define yynerrs         vmd_nerrs
#define yylloc          vmd_lloc

/* Copy the first part of user declarations.  */

/* Line 268 of yacc.c  */
#line 9 "vmd.y"

    #include "vmd.hxx"
    using namespace desres::msys::atomsel;


/* Line 268 of yacc.c  */
#line 85 "parser.cxx"

/* Enabling traces.  */
#ifndef YYDEBUG
# define YYDEBUG 0
#endif

/* Enabling verbose error messages.  */
#ifdef YYERROR_VERBOSE
# undef YYERROR_VERBOSE
# define YYERROR_VERBOSE 1
#else
# define YYERROR_VERBOSE 1
#endif

/* Enabling the token table.  */
#ifndef YYTOKEN_TABLE
# define YYTOKEN_TABLE 0
#endif


/* Tokens.  */
#ifndef YYTOKENTYPE
# define YYTOKENTYPE
   /* Put the tokens into the symbol table, so that GDB and other debuggers
      know about them.  */
   enum yytokentype {
     CMP = 258,
     INTKEY = 259,
     FLTKEY = 260,
     STRKEY = 261,
     WITHIN = 262,
     EXWITHIN = 263,
     PBWITHIN = 264,
     WITHINBONDS = 265,
     NEAREST = 266,
     SAME = 267,
     SINGLE = 268,
     IVAL = 269,
     FVAL = 270,
     SVAL = 271,
     FUNC = 272,
     REGEX = 273,
     MACRO = 274,
     STRFCTN = 275,
     AS = 276,
     OF = 277,
     TO = 278,
     ERR = 279,
     OR = 280,
     AND = 281,
     SUB = 282,
     ADD = 283,
     MOD = 284,
     DIV = 285,
     MUL = 286,
     EXP = 287,
     NOT = 288,
     nonassoc = 289
   };
#endif



#if ! defined YYSTYPE && ! defined YYSTYPE_IS_DECLARED
typedef union YYSTYPE
{

/* Line 301 of yacc.c  */
#line 14 "vmd.y"

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



/* Line 301 of yacc.c  */
#line 174 "parser.cxx"
} YYSTYPE;
# define YYSTYPE_IS_TRIVIAL 1
# define yystype YYSTYPE /* obsolescent; will be withdrawn */
# define YYSTYPE_IS_DECLARED 1
#endif

#if ! defined YYLTYPE && ! defined YYLTYPE_IS_DECLARED
typedef struct YYLTYPE
{
  int first_line;
  int first_column;
  int last_line;
  int last_column;
} YYLTYPE;
# define yyltype YYLTYPE /* obsolescent; will be withdrawn */
# define YYLTYPE_IS_DECLARED 1
# define YYLTYPE_IS_TRIVIAL 1
#endif


/* Copy the second part of user declarations.  */

/* Line 343 of yacc.c  */
#line 69 "vmd.y"

    int vmd_lex(YYSTYPE* lvalp, YYLTYPE* llocp, void* scanner);
    void vmd_error(YYLTYPE* locp, VMD* context, const char* err) {
      if (context->error.empty()) context->error = err;
    }

    #define scanner context->scanner


/* Line 343 of yacc.c  */
#line 209 "parser.cxx"

#ifdef short
# undef short
#endif

#ifdef YYTYPE_UINT8
typedef YYTYPE_UINT8 yytype_uint8;
#else
typedef unsigned char yytype_uint8;
#endif

#ifdef YYTYPE_INT8
typedef YYTYPE_INT8 yytype_int8;
#elif (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
typedef signed char yytype_int8;
#else
typedef short int yytype_int8;
#endif

#ifdef YYTYPE_UINT16
typedef YYTYPE_UINT16 yytype_uint16;
#else
typedef unsigned short int yytype_uint16;
#endif

#ifdef YYTYPE_INT16
typedef YYTYPE_INT16 yytype_int16;
#else
typedef short int yytype_int16;
#endif

#ifndef YYSIZE_T
# ifdef __SIZE_TYPE__
#  define YYSIZE_T __SIZE_TYPE__
# elif defined size_t
#  define YYSIZE_T size_t
# elif ! defined YYSIZE_T && (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
#  include <stddef.h> /* INFRINGES ON USER NAME SPACE */
#  define YYSIZE_T size_t
# else
#  define YYSIZE_T unsigned int
# endif
#endif

#define YYSIZE_MAXIMUM ((YYSIZE_T) -1)

#ifndef YY_
# if defined YYENABLE_NLS && YYENABLE_NLS
#  if ENABLE_NLS
#   include <libintl.h> /* INFRINGES ON USER NAME SPACE */
#   define YY_(msgid) dgettext ("bison-runtime", msgid)
#  endif
# endif
# ifndef YY_
#  define YY_(msgid) msgid
# endif
#endif

/* Suppress unused-variable warnings by "using" E.  */
#if ! defined lint || defined __GNUC__
# define YYUSE(e) ((void) (e))
#else
# define YYUSE(e) /* empty */
#endif

/* Identity function, used to suppress warnings about constant conditions.  */
#ifndef lint
# define YYID(n) (n)
#else
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static int
YYID (int yyi)
#else
static int
YYID (yyi)
    int yyi;
#endif
{
  return yyi;
}
#endif

#if ! defined yyoverflow || YYERROR_VERBOSE

/* The parser invokes alloca or malloc; define the necessary symbols.  */

# ifdef YYSTACK_USE_ALLOCA
#  if YYSTACK_USE_ALLOCA
#   ifdef __GNUC__
#    define YYSTACK_ALLOC __builtin_alloca
#   elif defined __BUILTIN_VA_ARG_INCR
#    include <alloca.h> /* INFRINGES ON USER NAME SPACE */
#   elif defined _AIX
#    define YYSTACK_ALLOC __alloca
#   elif defined _MSC_VER
#    include <malloc.h> /* INFRINGES ON USER NAME SPACE */
#    define alloca _alloca
#   else
#    define YYSTACK_ALLOC alloca
#    if ! defined _ALLOCA_H && ! defined EXIT_SUCCESS && (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
#     include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
#     ifndef EXIT_SUCCESS
#      define EXIT_SUCCESS 0
#     endif
#    endif
#   endif
#  endif
# endif

# ifdef YYSTACK_ALLOC
   /* Pacify GCC's `empty if-body' warning.  */
#  define YYSTACK_FREE(Ptr) do { /* empty */; } while (YYID (0))
#  ifndef YYSTACK_ALLOC_MAXIMUM
    /* The OS might guarantee only one guard page at the bottom of the stack,
       and a page size can be as small as 4096 bytes.  So we cannot safely
       invoke alloca (N) if N exceeds 4096.  Use a slightly smaller number
       to allow for a few compiler-allocated temporary stack slots.  */
#   define YYSTACK_ALLOC_MAXIMUM 4032 /* reasonable circa 2006 */
#  endif
# else
#  define YYSTACK_ALLOC YYMALLOC
#  define YYSTACK_FREE YYFREE
#  ifndef YYSTACK_ALLOC_MAXIMUM
#   define YYSTACK_ALLOC_MAXIMUM YYSIZE_MAXIMUM
#  endif
#  if (defined __cplusplus && ! defined EXIT_SUCCESS \
       && ! ((defined YYMALLOC || defined malloc) \
	     && (defined YYFREE || defined free)))
#   include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
#   ifndef EXIT_SUCCESS
#    define EXIT_SUCCESS 0
#   endif
#  endif
#  ifndef YYMALLOC
#   define YYMALLOC malloc
#   if ! defined malloc && ! defined EXIT_SUCCESS && (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
void *malloc (YYSIZE_T); /* INFRINGES ON USER NAME SPACE */
#   endif
#  endif
#  ifndef YYFREE
#   define YYFREE free
#   if ! defined free && ! defined EXIT_SUCCESS && (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
void free (void *); /* INFRINGES ON USER NAME SPACE */
#   endif
#  endif
# endif
#endif /* ! defined yyoverflow || YYERROR_VERBOSE */


#if (! defined yyoverflow \
     && (! defined __cplusplus \
	 || (defined YYLTYPE_IS_TRIVIAL && YYLTYPE_IS_TRIVIAL \
	     && defined YYSTYPE_IS_TRIVIAL && YYSTYPE_IS_TRIVIAL)))

/* A type that is properly aligned for any stack member.  */
union yyalloc
{
  yytype_int16 yyss_alloc;
  YYSTYPE yyvs_alloc;
  YYLTYPE yyls_alloc;
};

/* The size of the maximum gap between one aligned stack and the next.  */
# define YYSTACK_GAP_MAXIMUM (sizeof (union yyalloc) - 1)

/* The size of an array large to enough to hold all stacks, each with
   N elements.  */
# define YYSTACK_BYTES(N) \
     ((N) * (sizeof (yytype_int16) + sizeof (YYSTYPE) + sizeof (YYLTYPE)) \
      + 2 * YYSTACK_GAP_MAXIMUM)

# define YYCOPY_NEEDED 1

/* Relocate STACK from its old location to the new one.  The
   local variables YYSIZE and YYSTACKSIZE give the old and new number of
   elements in the stack, and YYPTR gives the new location of the
   stack.  Advance YYPTR to a properly aligned location for the next
   stack.  */
# define YYSTACK_RELOCATE(Stack_alloc, Stack)				\
    do									\
      {									\
	YYSIZE_T yynewbytes;						\
	YYCOPY (&yyptr->Stack_alloc, Stack, yysize);			\
	Stack = &yyptr->Stack_alloc;					\
	yynewbytes = yystacksize * sizeof (*Stack) + YYSTACK_GAP_MAXIMUM; \
	yyptr += yynewbytes / sizeof (*yyptr);				\
      }									\
    while (YYID (0))

#endif

#if defined YYCOPY_NEEDED && YYCOPY_NEEDED
/* Copy COUNT objects from FROM to TO.  The source and destination do
   not overlap.  */
# ifndef YYCOPY
#  if defined __GNUC__ && 1 < __GNUC__
#   define YYCOPY(To, From, Count) \
      __builtin_memcpy (To, From, (Count) * sizeof (*(From)))
#  else
#   define YYCOPY(To, From, Count)		\
      do					\
	{					\
	  YYSIZE_T yyi;				\
	  for (yyi = 0; yyi < (Count); yyi++)	\
	    (To)[yyi] = (From)[yyi];		\
	}					\
      while (YYID (0))
#  endif
# endif
#endif /* !YYCOPY_NEEDED */

/* YYFINAL -- State number of the termination state.  */
#define YYFINAL  53
/* YYLAST -- Last index in YYTABLE.  */
#define YYLAST   148

/* YYNTOKENS -- Number of terminals.  */
#define YYNTOKENS  37
/* YYNNTS -- Number of nonterminals.  */
#define YYNNTS  11
/* YYNRULES -- Number of rules.  */
#define YYNRULES  55
/* YYNRULES -- Number of states.  */
#define YYNSTATES  103

/* YYTRANSLATE(YYLEX) -- Bison symbol number corresponding to YYLEX.  */
#define YYUNDEFTOK  2
#define YYMAXUTOK   289

#define YYTRANSLATE(YYX)						\
  ((unsigned int) (YYX) <= YYMAXUTOK ? yytranslate[YYX] : YYUNDEFTOK)

/* YYTRANSLATE[YYLEX] -- Bison symbol number corresponding to YYLEX.  */
static const yytype_uint8 yytranslate[] =
{
       0,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
      35,    36,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     1,     2,     3,     4,
       5,     6,     7,     8,     9,    10,    11,    12,    13,    14,
      15,    16,    17,    18,    19,    20,    21,    22,    23,    24,
      25,    26,    27,    28,    29,    30,    31,    32,    33,    34
};

#if YYDEBUG
/* YYPRHS[YYN] -- Index of the first RHS symbol of rule number YYN in
   YYRHS.  */
static const yytype_uint8 yyprhs[] =
{
       0,     0,     3,     5,     7,    11,    14,    16,    18,    22,
      26,    28,    33,    38,    43,    48,    53,    58,    62,    64,
      66,    68,    71,    74,    77,    80,    82,    84,    86,    90,
      93,    98,   100,   104,   107,   112,   114,   116,   118,   121,
     124,   127,   129,   132,   134,   136,   138,   140,   143,   147,
     152,   156,   160,   164,   168,   172
};

/* YYRHS -- A `-1'-separated list of the rules' RHS.  */
static const yytype_int8 yyrhs[] =
{
      38,     0,    -1,    39,    -1,     1,    -1,    35,    39,    36,
      -1,    33,    39,    -1,    13,    -1,    19,    -1,    39,    26,
      39,    -1,    39,    25,    39,    -1,    41,    -1,     7,    42,
      22,    39,    -1,     8,    42,    22,    39,    -1,     9,    42,
      22,    39,    -1,    10,    14,    22,    39,    -1,    12,    40,
      21,    39,    -1,    11,    14,    23,    39,    -1,    47,     3,
      47,    -1,     4,    -1,     5,    -1,     6,    -1,     4,    43,
      -1,     5,    44,    -1,     6,    45,    -1,    20,    46,    -1,
      14,    -1,    15,    -1,    14,    -1,    14,    23,    14,    -1,
      43,    14,    -1,    43,    14,    23,    14,    -1,    15,    -1,
      15,    23,    15,    -1,    44,    15,    -1,    44,    15,    23,
      15,    -1,    16,    -1,    18,    -1,    14,    -1,    45,    16,
      -1,    45,    18,    -1,    45,    14,    -1,    16,    -1,    46,
      16,    -1,    15,    -1,    14,    -1,     5,    -1,     4,    -1,
      27,    47,    -1,    35,    47,    36,    -1,    17,    35,    47,
      36,    -1,    47,    28,    47,    -1,    47,    27,    47,    -1,
      47,    31,    47,    -1,    47,    30,    47,    -1,    47,    29,
      47,    -1,    47,    32,    47,    -1
};

/* YYRLINE[YYN] -- source line where rule number YYN was defined.  */
static const yytype_uint8 yyrline[] =
{
       0,    81,    81,    82,    86,    87,    88,    89,    90,    91,
      92,    93,    94,    95,    96,    97,    98,    99,   103,   104,
     105,   109,   110,   111,   112,   116,   117,   121,   122,   123,
     124,   128,   129,   130,   131,   135,   136,   137,   140,   141,
     142,   148,   149,   153,   154,   155,   156,   157,   158,   159,
     160,   161,   162,   163,   164,   165
};
#endif

#if YYDEBUG || YYERROR_VERBOSE || YYTOKEN_TABLE
/* YYTNAME[SYMBOL-NUM] -- String name of the symbol SYMBOL-NUM.
   First, the terminals, then, starting at YYNTOKENS, nonterminals.  */
static const char *const yytname[] =
{
  "$end", "error", "$undefined", "CMP", "INTKEY", "FLTKEY", "STRKEY",
  "WITHIN", "EXWITHIN", "PBWITHIN", "WITHINBONDS", "NEAREST", "SAME",
  "SINGLE", "IVAL", "FVAL", "SVAL", "FUNC", "REGEX", "MACRO", "STRFCTN",
  "AS", "OF", "TO", "ERR", "OR", "AND", "SUB", "ADD", "MOD", "DIV", "MUL",
  "EXP", "NOT", "nonassoc", "'('", "')'", "$accept", "start", "selection",
  "key", "key_value_list", "number", "intlist", "fltlist", "strlist",
  "stringlist", "expression", 0
};
#endif

# ifdef YYPRINT
/* YYTOKNUM[YYLEX-NUM] -- Internal token number corresponding to
   token YYLEX-NUM.  */
static const yytype_uint16 yytoknum[] =
{
       0,   256,   257,   258,   259,   260,   261,   262,   263,   264,
     265,   266,   267,   268,   269,   270,   271,   272,   273,   274,
     275,   276,   277,   278,   279,   280,   281,   282,   283,   284,
     285,   286,   287,   288,   289,    40,    41
};
# endif

/* YYR1[YYN] -- Symbol number of symbol that rule YYN derives.  */
static const yytype_uint8 yyr1[] =
{
       0,    37,    38,    38,    39,    39,    39,    39,    39,    39,
      39,    39,    39,    39,    39,    39,    39,    39,    40,    40,
      40,    41,    41,    41,    41,    42,    42,    43,    43,    43,
      43,    44,    44,    44,    44,    45,    45,    45,    45,    45,
      45,    46,    46,    47,    47,    47,    47,    47,    47,    47,
      47,    47,    47,    47,    47,    47
};

/* YYR2[YYN] -- Number of symbols composing right hand side of rule YYN.  */
static const yytype_uint8 yyr2[] =
{
       0,     2,     1,     1,     3,     2,     1,     1,     3,     3,
       1,     4,     4,     4,     4,     4,     4,     3,     1,     1,
       1,     2,     2,     2,     2,     1,     1,     1,     3,     2,
       4,     1,     3,     2,     4,     1,     1,     1,     2,     2,
       2,     1,     2,     1,     1,     1,     1,     2,     3,     4,
       3,     3,     3,     3,     3,     3
};

/* YYDEFACT[STATE-NAME] -- Default reduction number in state STATE-NUM.
   Performed when YYTABLE doesn't specify something else to do.  Zero
   means the default is an error.  */
static const yytype_uint8 yydefact[] =
{
       0,     3,    46,    45,     0,     0,     0,     0,     0,     0,
       0,     6,    44,    43,     0,     7,     0,     0,     0,     0,
       0,     2,    10,     0,    27,    21,    31,    22,    37,    35,
      36,    23,    25,    26,     0,     0,     0,     0,     0,    18,
      19,    20,     0,     0,    41,    24,    46,    45,     0,    47,
       5,     0,     0,     1,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,    29,     0,    33,    40,    38,    39,
       0,     0,     0,     0,     0,     0,     0,    42,     0,     4,
      48,     9,     8,    17,    51,    50,    54,    53,    52,    55,
      28,     0,    32,     0,    11,    12,    13,    14,    16,    15,
      49,    30,    34
};

/* YYDEFGOTO[NTERM-NUM].  */
static const yytype_int8 yydefgoto[] =
{
      -1,    20,    21,    42,    22,    34,    25,    27,    31,    45,
      23
};

/* YYPACT[STATE-NUM] -- Index in YYTABLE of the portion describing
   STATE-NUM.  */
#define YYPACT_NINF -19
static const yytype_int8 yypact[] =
{
       5,   -19,    -6,     6,   -11,    73,    73,    73,     9,    12,
      60,   -19,   -19,   -19,    -5,   -19,    23,    81,    64,    64,
      35,    77,   -19,    90,    57,    68,    69,    79,   -19,   -19,
     -19,    13,   -19,   -19,    82,    84,    85,    87,    89,   -19,
     -19,   -19,    92,    81,   -19,    94,   -19,   -19,    81,    19,
     -19,    75,    31,   -19,    64,    64,    81,    81,    81,    81,
      81,    81,    81,    91,   107,    99,   108,   -19,   -19,   -19,
      64,    64,    64,    64,    64,    64,    96,   -19,   106,   -19,
     -19,   103,   -19,   116,    19,    19,   109,   109,   109,   -19,
     -19,   101,   -19,   124,    77,    77,    77,    77,    77,    77,
     -19,   -19,   -19
};

/* YYPGOTO[NTERM-NUM].  */
static const yytype_int8 yypgoto[] =
{
     -19,   -19,   -18,   -19,   -19,    83,   -19,   -19,   -19,   -19,
     -15
};

/* YYTABLE[YYPACT[STATE-NUM]].  What to do in state STATE-NUM.  If
   positive, shift that token.  If negative, reduce the rule which
   number is the opposite.  If YYTABLE_NINF, syntax error.  */
#define YYTABLE_NINF -1
static const yytype_uint8 yytable[] =
{
      50,    51,    49,    28,    52,    29,     1,    30,    24,     2,
       3,     4,     5,     6,     7,     8,     9,    10,    11,    12,
      13,    26,    14,    37,    15,    16,    38,    67,    76,    68,
      43,    69,    17,    78,    56,    53,    81,    82,    18,    44,
      19,    83,    84,    85,    86,    87,    88,    89,    59,    60,
      61,    62,    94,    95,    96,    97,    98,    99,    57,    58,
      59,    60,    61,    62,    39,    40,    41,    80,     2,     3,
       4,     5,     6,     7,     8,     9,    10,    11,    12,    13,
      63,    14,    64,    15,    16,    46,    47,    32,    33,    35,
      36,    17,    65,    56,    66,    12,    13,    18,    14,    19,
      54,    55,    54,    55,    70,    90,    71,    72,    17,    73,
      77,    79,    74,    75,    92,   101,    48,    57,    58,    59,
      60,    61,    62,    57,    58,    59,    60,    61,    62,    55,
      91,    93,   100,    57,    58,    59,    60,    61,    62,   102,
       0,    62,    80,    57,    58,    59,    60,    61,    62
};

#define yypact_value_is_default(yystate) \
  ((yystate) == (-19))

#define yytable_value_is_error(yytable_value) \
  YYID (0)

static const yytype_int8 yycheck[] =
{
      18,    19,    17,    14,    19,    16,     1,    18,    14,     4,
       5,     6,     7,     8,     9,    10,    11,    12,    13,    14,
      15,    15,    17,    14,    19,    20,    14,    14,    43,    16,
      35,    18,    27,    48,     3,     0,    54,    55,    33,    16,
      35,    56,    57,    58,    59,    60,    61,    62,    29,    30,
      31,    32,    70,    71,    72,    73,    74,    75,    27,    28,
      29,    30,    31,    32,     4,     5,     6,    36,     4,     5,
       6,     7,     8,     9,    10,    11,    12,    13,    14,    15,
      23,    17,    14,    19,    20,     4,     5,    14,    15,     6,
       7,    27,    23,     3,    15,    14,    15,    33,    17,    35,
      25,    26,    25,    26,    22,    14,    22,    22,    27,    22,
      16,    36,    23,    21,    15,    14,    35,    27,    28,    29,
      30,    31,    32,    27,    28,    29,    30,    31,    32,    26,
      23,    23,    36,    27,    28,    29,    30,    31,    32,    15,
      -1,    32,    36,    27,    28,    29,    30,    31,    32
};

/* YYSTOS[STATE-NUM] -- The (internal number of the) accessing
   symbol of state STATE-NUM.  */
static const yytype_uint8 yystos[] =
{
       0,     1,     4,     5,     6,     7,     8,     9,    10,    11,
      12,    13,    14,    15,    17,    19,    20,    27,    33,    35,
      38,    39,    41,    47,    14,    43,    15,    44,    14,    16,
      18,    45,    14,    15,    42,    42,    42,    14,    14,     4,
       5,     6,    40,    35,    16,    46,     4,     5,    35,    47,
      39,    39,    47,     0,    25,    26,     3,    27,    28,    29,
      30,    31,    32,    23,    14,    23,    15,    14,    16,    18,
      22,    22,    22,    22,    23,    21,    47,    16,    47,    36,
      36,    39,    39,    47,    47,    47,    47,    47,    47,    47,
      14,    23,    15,    23,    39,    39,    39,    39,    39,    39,
      36,    14,    15
};

#define yyerrok		(yyerrstatus = 0)
#define yyclearin	(yychar = YYEMPTY)
#define YYEMPTY		(-2)
#define YYEOF		0

#define YYACCEPT	goto yyacceptlab
#define YYABORT		goto yyabortlab
#define YYERROR		goto yyerrorlab


/* Like YYERROR except do call yyerror.  This remains here temporarily
   to ease the transition to the new meaning of YYERROR, for GCC.
   Once GCC version 2 has supplanted version 1, this can go.  However,
   YYFAIL appears to be in use.  Nevertheless, it is formally deprecated
   in Bison 2.4.2's NEWS entry, where a plan to phase it out is
   discussed.  */

#define YYFAIL		goto yyerrlab
#if defined YYFAIL
  /* This is here to suppress warnings from the GCC cpp's
     -Wunused-macros.  Normally we don't worry about that warning, but
     some users do, and we want to make it easy for users to remove
     YYFAIL uses, which will produce warnings from Bison 2.5.  */
#endif

#define YYRECOVERING()  (!!yyerrstatus)

#define YYBACKUP(Token, Value)					\
do								\
  if (yychar == YYEMPTY && yylen == 1)				\
    {								\
      yychar = (Token);						\
      yylval = (Value);						\
      YYPOPSTACK (1);						\
      goto yybackup;						\
    }								\
  else								\
    {								\
      yyerror (&yylloc, context, YY_("syntax error: cannot back up")); \
      YYERROR;							\
    }								\
while (YYID (0))


#define YYTERROR	1
#define YYERRCODE	256


/* YYLLOC_DEFAULT -- Set CURRENT to span from RHS[1] to RHS[N].
   If N is 0, then set CURRENT to the empty location which ends
   the previous symbol: RHS[0] (always defined).  */

#define YYRHSLOC(Rhs, K) ((Rhs)[K])
#ifndef YYLLOC_DEFAULT
# define YYLLOC_DEFAULT(Current, Rhs, N)				\
    do									\
      if (YYID (N))                                                    \
	{								\
	  (Current).first_line   = YYRHSLOC (Rhs, 1).first_line;	\
	  (Current).first_column = YYRHSLOC (Rhs, 1).first_column;	\
	  (Current).last_line    = YYRHSLOC (Rhs, N).last_line;		\
	  (Current).last_column  = YYRHSLOC (Rhs, N).last_column;	\
	}								\
      else								\
	{								\
	  (Current).first_line   = (Current).last_line   =		\
	    YYRHSLOC (Rhs, 0).last_line;				\
	  (Current).first_column = (Current).last_column =		\
	    YYRHSLOC (Rhs, 0).last_column;				\
	}								\
    while (YYID (0))
#endif


/* YY_LOCATION_PRINT -- Print the location on the stream.
   This macro was not mandated originally: define only if we know
   we won't break user code: when these are the locations we know.  */

#ifndef YY_LOCATION_PRINT
# if defined YYLTYPE_IS_TRIVIAL && YYLTYPE_IS_TRIVIAL
#  define YY_LOCATION_PRINT(File, Loc)			\
     fprintf (File, "%d.%d-%d.%d",			\
	      (Loc).first_line, (Loc).first_column,	\
	      (Loc).last_line,  (Loc).last_column)
# else
#  define YY_LOCATION_PRINT(File, Loc) ((void) 0)
# endif
#endif


/* YYLEX -- calling `yylex' with the right arguments.  */

#ifdef YYLEX_PARAM
# define YYLEX yylex (&yylval, &yylloc, YYLEX_PARAM)
#else
# define YYLEX yylex (&yylval, &yylloc, scanner)
#endif

/* Enable debugging if requested.  */
#if YYDEBUG

# ifndef YYFPRINTF
#  include <stdio.h> /* INFRINGES ON USER NAME SPACE */
#  define YYFPRINTF fprintf
# endif

# define YYDPRINTF(Args)			\
do {						\
  if (yydebug)					\
    YYFPRINTF Args;				\
} while (YYID (0))

# define YY_SYMBOL_PRINT(Title, Type, Value, Location)			  \
do {									  \
  if (yydebug)								  \
    {									  \
      YYFPRINTF (stderr, "%s ", Title);					  \
      yy_symbol_print (stderr,						  \
		  Type, Value, Location, context); \
      YYFPRINTF (stderr, "\n");						  \
    }									  \
} while (YYID (0))


/*--------------------------------.
| Print this symbol on YYOUTPUT.  |
`--------------------------------*/

/*ARGSUSED*/
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yy_symbol_value_print (FILE *yyoutput, int yytype, YYSTYPE const * const yyvaluep, YYLTYPE const * const yylocationp, VMD* context)
#else
static void
yy_symbol_value_print (yyoutput, yytype, yyvaluep, yylocationp, context)
    FILE *yyoutput;
    int yytype;
    YYSTYPE const * const yyvaluep;
    YYLTYPE const * const yylocationp;
    VMD* context;
#endif
{
  if (!yyvaluep)
    return;
  YYUSE (yylocationp);
  YYUSE (context);
# ifdef YYPRINT
  if (yytype < YYNTOKENS)
    YYPRINT (yyoutput, yytoknum[yytype], *yyvaluep);
# else
  YYUSE (yyoutput);
# endif
  switch (yytype)
    {
      default:
	break;
    }
}


/*--------------------------------.
| Print this symbol on YYOUTPUT.  |
`--------------------------------*/

#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yy_symbol_print (FILE *yyoutput, int yytype, YYSTYPE const * const yyvaluep, YYLTYPE const * const yylocationp, VMD* context)
#else
static void
yy_symbol_print (yyoutput, yytype, yyvaluep, yylocationp, context)
    FILE *yyoutput;
    int yytype;
    YYSTYPE const * const yyvaluep;
    YYLTYPE const * const yylocationp;
    VMD* context;
#endif
{
  if (yytype < YYNTOKENS)
    YYFPRINTF (yyoutput, "token %s (", yytname[yytype]);
  else
    YYFPRINTF (yyoutput, "nterm %s (", yytname[yytype]);

  YY_LOCATION_PRINT (yyoutput, *yylocationp);
  YYFPRINTF (yyoutput, ": ");
  yy_symbol_value_print (yyoutput, yytype, yyvaluep, yylocationp, context);
  YYFPRINTF (yyoutput, ")");
}

/*------------------------------------------------------------------.
| yy_stack_print -- Print the state stack from its BOTTOM up to its |
| TOP (included).                                                   |
`------------------------------------------------------------------*/

#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yy_stack_print (yytype_int16 *yybottom, yytype_int16 *yytop)
#else
static void
yy_stack_print (yybottom, yytop)
    yytype_int16 *yybottom;
    yytype_int16 *yytop;
#endif
{
  YYFPRINTF (stderr, "Stack now");
  for (; yybottom <= yytop; yybottom++)
    {
      int yybot = *yybottom;
      YYFPRINTF (stderr, " %d", yybot);
    }
  YYFPRINTF (stderr, "\n");
}

# define YY_STACK_PRINT(Bottom, Top)				\
do {								\
  if (yydebug)							\
    yy_stack_print ((Bottom), (Top));				\
} while (YYID (0))


/*------------------------------------------------.
| Report that the YYRULE is going to be reduced.  |
`------------------------------------------------*/

#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yy_reduce_print (YYSTYPE *yyvsp, YYLTYPE *yylsp, int yyrule, VMD* context)
#else
static void
yy_reduce_print (yyvsp, yylsp, yyrule, context)
    YYSTYPE *yyvsp;
    YYLTYPE *yylsp;
    int yyrule;
    VMD* context;
#endif
{
  int yynrhs = yyr2[yyrule];
  int yyi;
  unsigned long int yylno = yyrline[yyrule];
  YYFPRINTF (stderr, "Reducing stack by rule %d (line %lu):\n",
	     yyrule - 1, yylno);
  /* The symbols being reduced.  */
  for (yyi = 0; yyi < yynrhs; yyi++)
    {
      YYFPRINTF (stderr, "   $%d = ", yyi + 1);
      yy_symbol_print (stderr, yyrhs[yyprhs[yyrule] + yyi],
		       &(yyvsp[(yyi + 1) - (yynrhs)])
		       , &(yylsp[(yyi + 1) - (yynrhs)])		       , context);
      YYFPRINTF (stderr, "\n");
    }
}

# define YY_REDUCE_PRINT(Rule)		\
do {					\
  if (yydebug)				\
    yy_reduce_print (yyvsp, yylsp, Rule, context); \
} while (YYID (0))

/* Nonzero means print parse trace.  It is left uninitialized so that
   multiple parsers can coexist.  */
int yydebug;
#else /* !YYDEBUG */
# define YYDPRINTF(Args)
# define YY_SYMBOL_PRINT(Title, Type, Value, Location)
# define YY_STACK_PRINT(Bottom, Top)
# define YY_REDUCE_PRINT(Rule)
#endif /* !YYDEBUG */


/* YYINITDEPTH -- initial size of the parser's stacks.  */
#ifndef	YYINITDEPTH
# define YYINITDEPTH 200
#endif

/* YYMAXDEPTH -- maximum size the stacks can grow to (effective only
   if the built-in stack extension method is used).

   Do not make this value too large; the results are undefined if
   YYSTACK_ALLOC_MAXIMUM < YYSTACK_BYTES (YYMAXDEPTH)
   evaluated with infinite-precision integer arithmetic.  */

#ifndef YYMAXDEPTH
# define YYMAXDEPTH 10000
#endif


#if YYERROR_VERBOSE

# ifndef yystrlen
#  if defined __GLIBC__ && defined _STRING_H
#   define yystrlen strlen
#  else
/* Return the length of YYSTR.  */
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static YYSIZE_T
yystrlen (const char *yystr)
#else
static YYSIZE_T
yystrlen (yystr)
    const char *yystr;
#endif
{
  YYSIZE_T yylen;
  for (yylen = 0; yystr[yylen]; yylen++)
    continue;
  return yylen;
}
#  endif
# endif

# ifndef yystpcpy
#  if defined __GLIBC__ && defined _STRING_H && defined _GNU_SOURCE
#   define yystpcpy stpcpy
#  else
/* Copy YYSRC to YYDEST, returning the address of the terminating '\0' in
   YYDEST.  */
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static char *
yystpcpy (char *yydest, const char *yysrc)
#else
static char *
yystpcpy (yydest, yysrc)
    char *yydest;
    const char *yysrc;
#endif
{
  char *yyd = yydest;
  const char *yys = yysrc;

  while ((*yyd++ = *yys++) != '\0')
    continue;

  return yyd - 1;
}
#  endif
# endif

# ifndef yytnamerr
/* Copy to YYRES the contents of YYSTR after stripping away unnecessary
   quotes and backslashes, so that it's suitable for yyerror.  The
   heuristic is that double-quoting is unnecessary unless the string
   contains an apostrophe, a comma, or backslash (other than
   backslash-backslash).  YYSTR is taken from yytname.  If YYRES is
   null, do not copy; instead, return the length of what the result
   would have been.  */
static YYSIZE_T
yytnamerr (char *yyres, const char *yystr)
{
  if (*yystr == '"')
    {
      YYSIZE_T yyn = 0;
      char const *yyp = yystr;

      for (;;)
	switch (*++yyp)
	  {
	  case '\'':
	  case ',':
	    goto do_not_strip_quotes;

	  case '\\':
	    if (*++yyp != '\\')
	      goto do_not_strip_quotes;
	    /* Fall through.  */
	  default:
	    if (yyres)
	      yyres[yyn] = *yyp;
	    yyn++;
	    break;

	  case '"':
	    if (yyres)
	      yyres[yyn] = '\0';
	    return yyn;
	  }
    do_not_strip_quotes: ;
    }

  if (! yyres)
    return yystrlen (yystr);

  return yystpcpy (yyres, yystr) - yyres;
}
# endif

/* Copy into *YYMSG, which is of size *YYMSG_ALLOC, an error message
   about the unexpected token YYTOKEN for the state stack whose top is
   YYSSP.

   Return 0 if *YYMSG was successfully written.  Return 1 if *YYMSG is
   not large enough to hold the message.  In that case, also set
   *YYMSG_ALLOC to the required number of bytes.  Return 2 if the
   required number of bytes is too large to store.  */
static int
yysyntax_error (YYSIZE_T *yymsg_alloc, char **yymsg,
                yytype_int16 *yyssp, int yytoken)
{
  YYSIZE_T yysize0 = yytnamerr (0, yytname[yytoken]);
  YYSIZE_T yysize = yysize0;
  YYSIZE_T yysize1;
  enum { YYERROR_VERBOSE_ARGS_MAXIMUM = 5 };
  /* Internationalized format string. */
  const char *yyformat = 0;
  /* Arguments of yyformat. */
  char const *yyarg[YYERROR_VERBOSE_ARGS_MAXIMUM];
  /* Number of reported tokens (one for the "unexpected", one per
     "expected"). */
  int yycount = 0;

  /* There are many possibilities here to consider:
     - Assume YYFAIL is not used.  It's too flawed to consider.  See
       <http://lists.gnu.org/archive/html/bison-patches/2009-12/msg00024.html>
       for details.  YYERROR is fine as it does not invoke this
       function.
     - If this state is a consistent state with a default action, then
       the only way this function was invoked is if the default action
       is an error action.  In that case, don't check for expected
       tokens because there are none.
     - The only way there can be no lookahead present (in yychar) is if
       this state is a consistent state with a default action.  Thus,
       detecting the absence of a lookahead is sufficient to determine
       that there is no unexpected or expected token to report.  In that
       case, just report a simple "syntax error".
     - Don't assume there isn't a lookahead just because this state is a
       consistent state with a default action.  There might have been a
       previous inconsistent state, consistent state with a non-default
       action, or user semantic action that manipulated yychar.
     - Of course, the expected token list depends on states to have
       correct lookahead information, and it depends on the parser not
       to perform extra reductions after fetching a lookahead from the
       scanner and before detecting a syntax error.  Thus, state merging
       (from LALR or IELR) and default reductions corrupt the expected
       token list.  However, the list is correct for canonical LR with
       one exception: it will still contain any token that will not be
       accepted due to an error action in a later state.
  */
  if (yytoken != YYEMPTY)
    {
      int yyn = yypact[*yyssp];
      yyarg[yycount++] = yytname[yytoken];
      if (!yypact_value_is_default (yyn))
        {
          /* Start YYX at -YYN if negative to avoid negative indexes in
             YYCHECK.  In other words, skip the first -YYN actions for
             this state because they are default actions.  */
          int yyxbegin = yyn < 0 ? -yyn : 0;
          /* Stay within bounds of both yycheck and yytname.  */
          int yychecklim = YYLAST - yyn + 1;
          int yyxend = yychecklim < YYNTOKENS ? yychecklim : YYNTOKENS;
          int yyx;

          for (yyx = yyxbegin; yyx < yyxend; ++yyx)
            if (yycheck[yyx + yyn] == yyx && yyx != YYTERROR
                && !yytable_value_is_error (yytable[yyx + yyn]))
              {
                if (yycount == YYERROR_VERBOSE_ARGS_MAXIMUM)
                  {
                    yycount = 1;
                    yysize = yysize0;
                    break;
                  }
                yyarg[yycount++] = yytname[yyx];
                yysize1 = yysize + yytnamerr (0, yytname[yyx]);
                if (! (yysize <= yysize1
                       && yysize1 <= YYSTACK_ALLOC_MAXIMUM))
                  return 2;
                yysize = yysize1;
              }
        }
    }

  switch (yycount)
    {
# define YYCASE_(N, S)                      \
      case N:                               \
        yyformat = S;                       \
      break
      YYCASE_(0, YY_("syntax error"));
      YYCASE_(1, YY_("syntax error, unexpected %s"));
      YYCASE_(2, YY_("syntax error, unexpected %s, expecting %s"));
      YYCASE_(3, YY_("syntax error, unexpected %s, expecting %s or %s"));
      YYCASE_(4, YY_("syntax error, unexpected %s, expecting %s or %s or %s"));
      YYCASE_(5, YY_("syntax error, unexpected %s, expecting %s or %s or %s or %s"));
# undef YYCASE_
    }

  yysize1 = yysize + yystrlen (yyformat);
  if (! (yysize <= yysize1 && yysize1 <= YYSTACK_ALLOC_MAXIMUM))
    return 2;
  yysize = yysize1;

  if (*yymsg_alloc < yysize)
    {
      *yymsg_alloc = 2 * yysize;
      if (! (yysize <= *yymsg_alloc
             && *yymsg_alloc <= YYSTACK_ALLOC_MAXIMUM))
        *yymsg_alloc = YYSTACK_ALLOC_MAXIMUM;
      return 1;
    }

  /* Avoid sprintf, as that infringes on the user's name space.
     Don't have undefined behavior even if the translation
     produced a string with the wrong number of "%s"s.  */
  {
    char *yyp = *yymsg;
    int yyi = 0;
    while ((*yyp = *yyformat) != '\0')
      if (*yyp == '%' && yyformat[1] == 's' && yyi < yycount)
        {
          yyp += yytnamerr (yyp, yyarg[yyi++]);
          yyformat += 2;
        }
      else
        {
          yyp++;
          yyformat++;
        }
  }
  return 0;
}
#endif /* YYERROR_VERBOSE */

/*-----------------------------------------------.
| Release the memory associated to this symbol.  |
`-----------------------------------------------*/

/*ARGSUSED*/
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yydestruct (const char *yymsg, int yytype, YYSTYPE *yyvaluep, YYLTYPE *yylocationp, VMD* context)
#else
static void
yydestruct (yymsg, yytype, yyvaluep, yylocationp, context)
    const char *yymsg;
    int yytype;
    YYSTYPE *yyvaluep;
    YYLTYPE *yylocationp;
    VMD* context;
#endif
{
  YYUSE (yyvaluep);
  YYUSE (yylocationp);
  YYUSE (context);

  if (!yymsg)
    yymsg = "Deleting";
  YY_SYMBOL_PRINT (yymsg, yytype, yyvaluep, yylocationp);

  switch (yytype)
    {
      case 16: /* "SVAL" */

/* Line 1391 of yacc.c  */
#line 66 "vmd.y"
	{ free((yyvaluep->sval)); };

/* Line 1391 of yacc.c  */
#line 1257 "parser.cxx"
	break;
      case 41: /* "key_value_list" */

/* Line 1391 of yacc.c  */
#line 67 "vmd.y"
	{delete (yyvaluep->pred); };

/* Line 1391 of yacc.c  */
#line 1266 "parser.cxx"
	break;
      case 43: /* "intlist" */

/* Line 1391 of yacc.c  */
#line 67 "vmd.y"
	{delete (yyvaluep->ilist); };

/* Line 1391 of yacc.c  */
#line 1275 "parser.cxx"
	break;
      case 44: /* "fltlist" */

/* Line 1391 of yacc.c  */
#line 67 "vmd.y"
	{delete (yyvaluep->flist); };

/* Line 1391 of yacc.c  */
#line 1284 "parser.cxx"
	break;
      case 45: /* "strlist" */

/* Line 1391 of yacc.c  */
#line 67 "vmd.y"
	{delete (yyvaluep->slist); };

/* Line 1391 of yacc.c  */
#line 1293 "parser.cxx"
	break;

      default:
	break;
    }
}


/* Prevent warnings from -Wmissing-prototypes.  */
#ifdef YYPARSE_PARAM
#if defined __STDC__ || defined __cplusplus
int yyparse (void *YYPARSE_PARAM);
#else
int yyparse ();
#endif
#else /* ! YYPARSE_PARAM */
#if defined __STDC__ || defined __cplusplus
int yyparse (VMD* context);
#else
int yyparse ();
#endif
#endif /* ! YYPARSE_PARAM */


/*----------.
| yyparse.  |
`----------*/

#ifdef YYPARSE_PARAM
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
int
yyparse (void *YYPARSE_PARAM)
#else
int
yyparse (YYPARSE_PARAM)
    void *YYPARSE_PARAM;
#endif
#else /* ! YYPARSE_PARAM */
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
int
yyparse (VMD* context)
#else
int
yyparse (context)
    VMD* context;
#endif
#endif
{
/* The lookahead symbol.  */
int yychar;

/* The semantic value of the lookahead symbol.  */
YYSTYPE yylval;

/* Location data for the lookahead symbol.  */
YYLTYPE yylloc;

    /* Number of syntax errors so far.  */
    int yynerrs;

    int yystate;
    /* Number of tokens to shift before error messages enabled.  */
    int yyerrstatus;

    /* The stacks and their tools:
       `yyss': related to states.
       `yyvs': related to semantic values.
       `yyls': related to locations.

       Refer to the stacks thru separate pointers, to allow yyoverflow
       to reallocate them elsewhere.  */

    /* The state stack.  */
    yytype_int16 yyssa[YYINITDEPTH];
    yytype_int16 *yyss;
    yytype_int16 *yyssp;

    /* The semantic value stack.  */
    YYSTYPE yyvsa[YYINITDEPTH];
    YYSTYPE *yyvs;
    YYSTYPE *yyvsp;

    /* The location stack.  */
    YYLTYPE yylsa[YYINITDEPTH];
    YYLTYPE *yyls;
    YYLTYPE *yylsp;

    /* The locations where the error started and ended.  */
    YYLTYPE yyerror_range[3];

    YYSIZE_T yystacksize;

  int yyn;
  int yyresult;
  /* Lookahead token as an internal (translated) token number.  */
  int yytoken;
  /* The variables used to return semantic value and location from the
     action routines.  */
  YYSTYPE yyval;
  YYLTYPE yyloc;

#if YYERROR_VERBOSE
  /* Buffer for error messages, and its allocated size.  */
  char yymsgbuf[128];
  char *yymsg = yymsgbuf;
  YYSIZE_T yymsg_alloc = sizeof yymsgbuf;
#endif

#define YYPOPSTACK(N)   (yyvsp -= (N), yyssp -= (N), yylsp -= (N))

  /* The number of symbols on the RHS of the reduced rule.
     Keep to zero when no symbol should be popped.  */
  int yylen = 0;

  yytoken = 0;
  yyss = yyssa;
  yyvs = yyvsa;
  yyls = yylsa;
  yystacksize = YYINITDEPTH;

  YYDPRINTF ((stderr, "Starting parse\n"));

  yystate = 0;
  yyerrstatus = 0;
  yynerrs = 0;
  yychar = YYEMPTY; /* Cause a token to be read.  */

  /* Initialize stack pointers.
     Waste one element of value and location stack
     so that they stay on the same level as the state stack.
     The wasted elements are never initialized.  */
  yyssp = yyss;
  yyvsp = yyvs;
  yylsp = yyls;

#if defined YYLTYPE_IS_TRIVIAL && YYLTYPE_IS_TRIVIAL
  /* Initialize the default location before parsing starts.  */
  yylloc.first_line   = yylloc.last_line   = 1;
  yylloc.first_column = yylloc.last_column = 1;
#endif

  goto yysetstate;

/*------------------------------------------------------------.
| yynewstate -- Push a new state, which is found in yystate.  |
`------------------------------------------------------------*/
 yynewstate:
  /* In all cases, when you get here, the value and location stacks
     have just been pushed.  So pushing a state here evens the stacks.  */
  yyssp++;

 yysetstate:
  *yyssp = yystate;

  if (yyss + yystacksize - 1 <= yyssp)
    {
      /* Get the current used size of the three stacks, in elements.  */
      YYSIZE_T yysize = yyssp - yyss + 1;

#ifdef yyoverflow
      {
	/* Give user a chance to reallocate the stack.  Use copies of
	   these so that the &'s don't force the real ones into
	   memory.  */
	YYSTYPE *yyvs1 = yyvs;
	yytype_int16 *yyss1 = yyss;
	YYLTYPE *yyls1 = yyls;

	/* Each stack pointer address is followed by the size of the
	   data in use in that stack, in bytes.  This used to be a
	   conditional around just the two extra args, but that might
	   be undefined if yyoverflow is a macro.  */
	yyoverflow (YY_("memory exhausted"),
		    &yyss1, yysize * sizeof (*yyssp),
		    &yyvs1, yysize * sizeof (*yyvsp),
		    &yyls1, yysize * sizeof (*yylsp),
		    &yystacksize);

	yyls = yyls1;
	yyss = yyss1;
	yyvs = yyvs1;
      }
#else /* no yyoverflow */
# ifndef YYSTACK_RELOCATE
      goto yyexhaustedlab;
# else
      /* Extend the stack our own way.  */
      if (YYMAXDEPTH <= yystacksize)
	goto yyexhaustedlab;
      yystacksize *= 2;
      if (YYMAXDEPTH < yystacksize)
	yystacksize = YYMAXDEPTH;

      {
	yytype_int16 *yyss1 = yyss;
	union yyalloc *yyptr =
	  (union yyalloc *) YYSTACK_ALLOC (YYSTACK_BYTES (yystacksize));
	if (! yyptr)
	  goto yyexhaustedlab;
	YYSTACK_RELOCATE (yyss_alloc, yyss);
	YYSTACK_RELOCATE (yyvs_alloc, yyvs);
	YYSTACK_RELOCATE (yyls_alloc, yyls);
#  undef YYSTACK_RELOCATE
	if (yyss1 != yyssa)
	  YYSTACK_FREE (yyss1);
      }
# endif
#endif /* no yyoverflow */

      yyssp = yyss + yysize - 1;
      yyvsp = yyvs + yysize - 1;
      yylsp = yyls + yysize - 1;

      YYDPRINTF ((stderr, "Stack size increased to %lu\n",
		  (unsigned long int) yystacksize));

      if (yyss + yystacksize - 1 <= yyssp)
	YYABORT;
    }

  YYDPRINTF ((stderr, "Entering state %d\n", yystate));

  if (yystate == YYFINAL)
    YYACCEPT;

  goto yybackup;

/*-----------.
| yybackup.  |
`-----------*/
yybackup:

  /* Do appropriate processing given the current state.  Read a
     lookahead token if we need one and don't already have one.  */

  /* First try to decide what to do without reference to lookahead token.  */
  yyn = yypact[yystate];
  if (yypact_value_is_default (yyn))
    goto yydefault;

  /* Not known => get a lookahead token if don't already have one.  */

  /* YYCHAR is either YYEMPTY or YYEOF or a valid lookahead symbol.  */
  if (yychar == YYEMPTY)
    {
      YYDPRINTF ((stderr, "Reading a token: "));
      yychar = YYLEX;
    }

  if (yychar <= YYEOF)
    {
      yychar = yytoken = YYEOF;
      YYDPRINTF ((stderr, "Now at end of input.\n"));
    }
  else
    {
      yytoken = YYTRANSLATE (yychar);
      YY_SYMBOL_PRINT ("Next token is", yytoken, &yylval, &yylloc);
    }

  /* If the proper action on seeing token YYTOKEN is to reduce or to
     detect an error, take that action.  */
  yyn += yytoken;
  if (yyn < 0 || YYLAST < yyn || yycheck[yyn] != yytoken)
    goto yydefault;
  yyn = yytable[yyn];
  if (yyn <= 0)
    {
      if (yytable_value_is_error (yyn))
        goto yyerrlab;
      yyn = -yyn;
      goto yyreduce;
    }

  /* Count tokens shifted since error; after three, turn off error
     status.  */
  if (yyerrstatus)
    yyerrstatus--;

  /* Shift the lookahead token.  */
  YY_SYMBOL_PRINT ("Shifting", yytoken, &yylval, &yylloc);

  /* Discard the shifted token.  */
  yychar = YYEMPTY;

  yystate = yyn;
  *++yyvsp = yylval;
  *++yylsp = yylloc;
  goto yynewstate;


/*-----------------------------------------------------------.
| yydefault -- do the default action for the current state.  |
`-----------------------------------------------------------*/
yydefault:
  yyn = yydefact[yystate];
  if (yyn == 0)
    goto yyerrlab;
  goto yyreduce;


/*-----------------------------.
| yyreduce -- Do a reduction.  |
`-----------------------------*/
yyreduce:
  /* yyn is the number of a rule to reduce with.  */
  yylen = yyr2[yyn];

  /* If YYLEN is nonzero, implement the default value of the action:
     `$$ = $1'.

     Otherwise, the following line sets YYVAL to garbage.
     This behavior is undocumented and Bison
     users should not rely upon it.  Assigning to YYVAL
     unconditionally makes the parser a bit smaller, and it avoids a
     GCC warning that YYVAL may be used uninitialized.  */
  yyval = yyvsp[1-yylen];

  /* Default location.  */
  YYLLOC_DEFAULT (yyloc, (yylsp - yylen), yylen);
  YY_REDUCE_PRINT (yyn);
  switch (yyn)
    {
        case 2:

/* Line 1821 of yacc.c  */
#line 81 "vmd.y"
    { context->result = (yyvsp[(1) - (1)].pred); }
    break;

  case 3:

/* Line 1821 of yacc.c  */
#line 82 "vmd.y"
    { context->result = NULL; }
    break;

  case 4:

/* Line 1821 of yacc.c  */
#line 86 "vmd.y"
    { (yyval.pred) = (yyvsp[(2) - (3)].pred);                       }
    break;

  case 5:

/* Line 1821 of yacc.c  */
#line 87 "vmd.y"
    { (yyval.pred) = context->make_not((yyvsp[(2) - (2)].pred));    }
    break;

  case 6:

/* Line 1821 of yacc.c  */
#line 88 "vmd.y"
    { (yyval.pred) = context->make_single((yyvsp[(1) - (1)].key)); }
    break;

  case 7:

/* Line 1821 of yacc.c  */
#line 89 "vmd.y"
    { (yyval.pred) = (yyvsp[(1) - (1)].pred); }
    break;

  case 8:

/* Line 1821 of yacc.c  */
#line 90 "vmd.y"
    { (yyval.pred) = context->make_and((yyvsp[(1) - (3)].pred),(yyvsp[(3) - (3)].pred)); }
    break;

  case 9:

/* Line 1821 of yacc.c  */
#line 91 "vmd.y"
    { (yyval.pred) = context->make_or((yyvsp[(1) - (3)].pred),(yyvsp[(3) - (3)].pred));  }
    break;

  case 10:

/* Line 1821 of yacc.c  */
#line 92 "vmd.y"
    { (yyval.pred) = (yyvsp[(1) - (1)].pred); }
    break;

  case 11:

/* Line 1821 of yacc.c  */
#line 93 "vmd.y"
    { (yyval.pred) = context->make_within((yyvsp[(2) - (4)].fval),(yyvsp[(4) - (4)].pred)); }
    break;

  case 12:

/* Line 1821 of yacc.c  */
#line 94 "vmd.y"
    { (yyval.pred) = context->make_exwithin((yyvsp[(2) - (4)].fval),(yyvsp[(4) - (4)].pred)); }
    break;

  case 13:

/* Line 1821 of yacc.c  */
#line 95 "vmd.y"
    { (yyval.pred) = context->make_pbwithin((yyvsp[(2) - (4)].fval),(yyvsp[(4) - (4)].pred)); }
    break;

  case 14:

/* Line 1821 of yacc.c  */
#line 96 "vmd.y"
    { (yyval.pred) = context->make_withinbonds((yyvsp[(2) - (4)].ival),(yyvsp[(4) - (4)].pred)); }
    break;

  case 15:

/* Line 1821 of yacc.c  */
#line 97 "vmd.y"
    { (yyval.pred)= context->make_same((yyvsp[(2) - (4)].key),(yyvsp[(4) - (4)].pred)); }
    break;

  case 16:

/* Line 1821 of yacc.c  */
#line 98 "vmd.y"
    { (yyval.pred) = context->make_nearest((yyvsp[(2) - (4)].ival),(yyvsp[(4) - (4)].pred)); }
    break;

  case 17:

/* Line 1821 of yacc.c  */
#line 99 "vmd.y"
    { (yyval.pred) = context->make_compare((yyvsp[(2) - (3)].ival),(yyvsp[(1) - (3)].expr),(yyvsp[(3) - (3)].expr)); }
    break;

  case 18:

/* Line 1821 of yacc.c  */
#line 103 "vmd.y"
    { (yyval.key) = (yyvsp[(1) - (1)].key); }
    break;

  case 19:

/* Line 1821 of yacc.c  */
#line 104 "vmd.y"
    { (yyval.key) = (yyvsp[(1) - (1)].key); }
    break;

  case 20:

/* Line 1821 of yacc.c  */
#line 105 "vmd.y"
    { (yyval.key) = (yyvsp[(1) - (1)].key); }
    break;

  case 21:

/* Line 1821 of yacc.c  */
#line 109 "vmd.y"
    { (yyval.pred) = context->make_key((yyvsp[(1) - (2)].key),(yyvsp[(2) - (2)].ilist)); }
    break;

  case 22:

/* Line 1821 of yacc.c  */
#line 110 "vmd.y"
    { (yyval.pred) = context->make_key((yyvsp[(1) - (2)].key),(yyvsp[(2) - (2)].flist)); }
    break;

  case 23:

/* Line 1821 of yacc.c  */
#line 111 "vmd.y"
    { (yyval.pred) = context->make_key((yyvsp[(1) - (2)].key),(yyvsp[(2) - (2)].slist)); }
    break;

  case 24:

/* Line 1821 of yacc.c  */
#line 112 "vmd.y"
    { (yyval.pred) = context->make_strfctn((yyvsp[(1) - (2)].strfctn),(yyvsp[(2) - (2)].slist)); }
    break;

  case 25:

/* Line 1821 of yacc.c  */
#line 116 "vmd.y"
    { (yyval.fval) = (yyvsp[(1) - (1)].ival); }
    break;

  case 26:

/* Line 1821 of yacc.c  */
#line 117 "vmd.y"
    { (yyval.fval) = (yyvsp[(1) - (1)].fval); }
    break;

  case 27:

/* Line 1821 of yacc.c  */
#line 121 "vmd.y"
    { (yyval.ilist) = new IntList((yyvsp[(1) - (1)].ival)); }
    break;

  case 28:

/* Line 1821 of yacc.c  */
#line 122 "vmd.y"
    { (yyval.ilist) = new IntList((yyvsp[(1) - (3)].ival),(yyvsp[(3) - (3)].ival)); }
    break;

  case 29:

/* Line 1821 of yacc.c  */
#line 123 "vmd.y"
    { (yyvsp[(1) - (2)].ilist)->add((yyvsp[(2) - (2)].ival)); (yyval.ilist)=(yyvsp[(1) - (2)].ilist); }
    break;

  case 30:

/* Line 1821 of yacc.c  */
#line 124 "vmd.y"
    { (yyvsp[(1) - (4)].ilist)->add((yyvsp[(2) - (4)].ival),(yyvsp[(4) - (4)].ival)); (yyval.ilist)=(yyvsp[(1) - (4)].ilist); }
    break;

  case 31:

/* Line 1821 of yacc.c  */
#line 128 "vmd.y"
    { (yyval.flist) = new FloatList((yyvsp[(1) - (1)].fval)); }
    break;

  case 32:

/* Line 1821 of yacc.c  */
#line 129 "vmd.y"
    { (yyval.flist) = new FloatList((yyvsp[(1) - (3)].fval),(yyvsp[(3) - (3)].fval)); }
    break;

  case 33:

/* Line 1821 of yacc.c  */
#line 130 "vmd.y"
    { (yyvsp[(1) - (2)].flist)->add((yyvsp[(2) - (2)].fval)); (yyval.flist)=(yyvsp[(1) - (2)].flist); }
    break;

  case 34:

/* Line 1821 of yacc.c  */
#line 131 "vmd.y"
    { (yyvsp[(1) - (4)].flist)->add((yyvsp[(2) - (4)].fval),(yyvsp[(4) - (4)].fval)); (yyval.flist)=(yyvsp[(1) - (4)].flist); }
    break;

  case 35:

/* Line 1821 of yacc.c  */
#line 135 "vmd.y"
    { (yyval.slist) = new StringList((yyvsp[(1) - (1)].sval)); free((yyvsp[(1) - (1)].sval)); }
    break;

  case 36:

/* Line 1821 of yacc.c  */
#line 136 "vmd.y"
    { (yyval.slist) = new StringList((yyvsp[(1) - (1)].regex)); }
    break;

  case 37:

/* Line 1821 of yacc.c  */
#line 137 "vmd.y"
    { char buf[32];
                                  sprintf(buf, "%d", (yyvsp[(1) - (1)].ival)); 
                                  (yyval.slist) = new StringList(buf); }
    break;

  case 38:

/* Line 1821 of yacc.c  */
#line 140 "vmd.y"
    { (yyvsp[(1) - (2)].slist)->add((yyvsp[(2) - (2)].sval)); free((yyvsp[(2) - (2)].sval)); (yyval.slist)=(yyvsp[(1) - (2)].slist);      }
    break;

  case 39:

/* Line 1821 of yacc.c  */
#line 141 "vmd.y"
    { (yyvsp[(1) - (2)].slist)->add((yyvsp[(2) - (2)].regex)); (yyval.slist)=(yyvsp[(1) - (2)].slist);      }
    break;

  case 40:

/* Line 1821 of yacc.c  */
#line 142 "vmd.y"
    { char buf[32];
                                  sprintf(buf, "%d", (yyvsp[(2) - (2)].ival)); 
                                  (yyvsp[(1) - (2)].slist)->add(buf); (yyval.slist)=(yyvsp[(1) - (2)].slist);     }
    break;

  case 41:

/* Line 1821 of yacc.c  */
#line 148 "vmd.y"
    { (yyval.slist) = new StringList((yyvsp[(1) - (1)].sval)); free((yyvsp[(1) - (1)].sval)); }
    break;

  case 42:

/* Line 1821 of yacc.c  */
#line 149 "vmd.y"
    { (yyvsp[(1) - (2)].slist)->add((yyvsp[(2) - (2)].sval)); free((yyvsp[(2) - (2)].sval)); (yyval.slist)=(yyvsp[(1) - (2)].slist); }
    break;

  case 43:

/* Line 1821 of yacc.c  */
#line 153 "vmd.y"
    { (yyval.expr) = context->make_exp((yyvsp[(1) - (1)].fval)); }
    break;

  case 44:

/* Line 1821 of yacc.c  */
#line 154 "vmd.y"
    { (yyval.expr) = context->make_exp((yyvsp[(1) - (1)].ival)); }
    break;

  case 45:

/* Line 1821 of yacc.c  */
#line 155 "vmd.y"
    { (yyval.expr) = context->make_exp((yyvsp[(1) - (1)].key)); }
    break;

  case 46:

/* Line 1821 of yacc.c  */
#line 156 "vmd.y"
    { (yyval.expr) = context->make_exp((yyvsp[(1) - (1)].key)); }
    break;

  case 47:

/* Line 1821 of yacc.c  */
#line 157 "vmd.y"
    { (yyval.expr) = context->make_unaexp("-", (yyvsp[(2) - (2)].expr)); }
    break;

  case 48:

/* Line 1821 of yacc.c  */
#line 158 "vmd.y"
    { (yyval.expr) = (yyvsp[(2) - (3)].expr); }
    break;

  case 49:

/* Line 1821 of yacc.c  */
#line 159 "vmd.y"
    { (yyval.expr) = context->make_unaexp((yyvsp[(1) - (4)].sval),(yyvsp[(3) - (4)].expr)); free((yyvsp[(1) - (4)].sval)); }
    break;

  case 50:

/* Line 1821 of yacc.c  */
#line 160 "vmd.y"
    { (yyval.expr) = context->make_binexp("+", (yyvsp[(1) - (3)].expr),(yyvsp[(3) - (3)].expr)); }
    break;

  case 51:

/* Line 1821 of yacc.c  */
#line 161 "vmd.y"
    { (yyval.expr) = context->make_binexp("-", (yyvsp[(1) - (3)].expr),(yyvsp[(3) - (3)].expr)); }
    break;

  case 52:

/* Line 1821 of yacc.c  */
#line 162 "vmd.y"
    { (yyval.expr) = context->make_binexp("*", (yyvsp[(1) - (3)].expr),(yyvsp[(3) - (3)].expr)); }
    break;

  case 53:

/* Line 1821 of yacc.c  */
#line 163 "vmd.y"
    { (yyval.expr) = context->make_binexp("/", (yyvsp[(1) - (3)].expr),(yyvsp[(3) - (3)].expr)); }
    break;

  case 54:

/* Line 1821 of yacc.c  */
#line 164 "vmd.y"
    { (yyval.expr) = context->make_binexp("%", (yyvsp[(1) - (3)].expr),(yyvsp[(3) - (3)].expr)); }
    break;

  case 55:

/* Line 1821 of yacc.c  */
#line 165 "vmd.y"
    { (yyval.expr) = context->make_binexp("^", (yyvsp[(1) - (3)].expr),(yyvsp[(3) - (3)].expr)); }
    break;



/* Line 1821 of yacc.c  */
#line 2004 "parser.cxx"
      default: break;
    }
  /* User semantic actions sometimes alter yychar, and that requires
     that yytoken be updated with the new translation.  We take the
     approach of translating immediately before every use of yytoken.
     One alternative is translating here after every semantic action,
     but that translation would be missed if the semantic action invokes
     YYABORT, YYACCEPT, or YYERROR immediately after altering yychar or
     if it invokes YYBACKUP.  In the case of YYABORT or YYACCEPT, an
     incorrect destructor might then be invoked immediately.  In the
     case of YYERROR or YYBACKUP, subsequent parser actions might lead
     to an incorrect destructor call or verbose syntax error message
     before the lookahead is translated.  */
  YY_SYMBOL_PRINT ("-> $$ =", yyr1[yyn], &yyval, &yyloc);

  YYPOPSTACK (yylen);
  yylen = 0;
  YY_STACK_PRINT (yyss, yyssp);

  *++yyvsp = yyval;
  *++yylsp = yyloc;

  /* Now `shift' the result of the reduction.  Determine what state
     that goes to, based on the state we popped back to and the rule
     number reduced by.  */

  yyn = yyr1[yyn];

  yystate = yypgoto[yyn - YYNTOKENS] + *yyssp;
  if (0 <= yystate && yystate <= YYLAST && yycheck[yystate] == *yyssp)
    yystate = yytable[yystate];
  else
    yystate = yydefgoto[yyn - YYNTOKENS];

  goto yynewstate;


/*------------------------------------.
| yyerrlab -- here on detecting error |
`------------------------------------*/
yyerrlab:
  /* Make sure we have latest lookahead translation.  See comments at
     user semantic actions for why this is necessary.  */
  yytoken = yychar == YYEMPTY ? YYEMPTY : YYTRANSLATE (yychar);

  /* If not already recovering from an error, report this error.  */
  if (!yyerrstatus)
    {
      ++yynerrs;
#if ! YYERROR_VERBOSE
      yyerror (&yylloc, context, YY_("syntax error"));
#else
# define YYSYNTAX_ERROR yysyntax_error (&yymsg_alloc, &yymsg, \
                                        yyssp, yytoken)
      {
        char const *yymsgp = YY_("syntax error");
        int yysyntax_error_status;
        yysyntax_error_status = YYSYNTAX_ERROR;
        if (yysyntax_error_status == 0)
          yymsgp = yymsg;
        else if (yysyntax_error_status == 1)
          {
            if (yymsg != yymsgbuf)
              YYSTACK_FREE (yymsg);
            yymsg = (char *) YYSTACK_ALLOC (yymsg_alloc);
            if (!yymsg)
              {
                yymsg = yymsgbuf;
                yymsg_alloc = sizeof yymsgbuf;
                yysyntax_error_status = 2;
              }
            else
              {
                yysyntax_error_status = YYSYNTAX_ERROR;
                yymsgp = yymsg;
              }
          }
        yyerror (&yylloc, context, yymsgp);
        if (yysyntax_error_status == 2)
          goto yyexhaustedlab;
      }
# undef YYSYNTAX_ERROR
#endif
    }

  yyerror_range[1] = yylloc;

  if (yyerrstatus == 3)
    {
      /* If just tried and failed to reuse lookahead token after an
	 error, discard it.  */

      if (yychar <= YYEOF)
	{
	  /* Return failure if at end of input.  */
	  if (yychar == YYEOF)
	    YYABORT;
	}
      else
	{
	  yydestruct ("Error: discarding",
		      yytoken, &yylval, &yylloc, context);
	  yychar = YYEMPTY;
	}
    }

  /* Else will try to reuse lookahead token after shifting the error
     token.  */
  goto yyerrlab1;


/*---------------------------------------------------.
| yyerrorlab -- error raised explicitly by YYERROR.  |
`---------------------------------------------------*/
yyerrorlab:

  /* Pacify compilers like GCC when the user code never invokes
     YYERROR and the label yyerrorlab therefore never appears in user
     code.  */
  if (/*CONSTCOND*/ 0)
     goto yyerrorlab;

  yyerror_range[1] = yylsp[1-yylen];
  /* Do not reclaim the symbols of the rule which action triggered
     this YYERROR.  */
  YYPOPSTACK (yylen);
  yylen = 0;
  YY_STACK_PRINT (yyss, yyssp);
  yystate = *yyssp;
  goto yyerrlab1;


/*-------------------------------------------------------------.
| yyerrlab1 -- common code for both syntax error and YYERROR.  |
`-------------------------------------------------------------*/
yyerrlab1:
  yyerrstatus = 3;	/* Each real token shifted decrements this.  */

  for (;;)
    {
      yyn = yypact[yystate];
      if (!yypact_value_is_default (yyn))
	{
	  yyn += YYTERROR;
	  if (0 <= yyn && yyn <= YYLAST && yycheck[yyn] == YYTERROR)
	    {
	      yyn = yytable[yyn];
	      if (0 < yyn)
		break;
	    }
	}

      /* Pop the current state because it cannot handle the error token.  */
      if (yyssp == yyss)
	YYABORT;

      yyerror_range[1] = *yylsp;
      yydestruct ("Error: popping",
		  yystos[yystate], yyvsp, yylsp, context);
      YYPOPSTACK (1);
      yystate = *yyssp;
      YY_STACK_PRINT (yyss, yyssp);
    }

  *++yyvsp = yylval;

  yyerror_range[2] = yylloc;
  /* Using YYLLOC is tempting, but would change the location of
     the lookahead.  YYLOC is available though.  */
  YYLLOC_DEFAULT (yyloc, yyerror_range, 2);
  *++yylsp = yyloc;

  /* Shift the error token.  */
  YY_SYMBOL_PRINT ("Shifting", yystos[yyn], yyvsp, yylsp);

  yystate = yyn;
  goto yynewstate;


/*-------------------------------------.
| yyacceptlab -- YYACCEPT comes here.  |
`-------------------------------------*/
yyacceptlab:
  yyresult = 0;
  goto yyreturn;

/*-----------------------------------.
| yyabortlab -- YYABORT comes here.  |
`-----------------------------------*/
yyabortlab:
  yyresult = 1;
  goto yyreturn;

#if !defined(yyoverflow) || YYERROR_VERBOSE
/*-------------------------------------------------.
| yyexhaustedlab -- memory exhaustion comes here.  |
`-------------------------------------------------*/
yyexhaustedlab:
  yyerror (&yylloc, context, YY_("memory exhausted"));
  yyresult = 2;
  /* Fall through.  */
#endif

yyreturn:
  if (yychar != YYEMPTY)
    {
      /* Make sure we have latest lookahead translation.  See comments at
         user semantic actions for why this is necessary.  */
      yytoken = YYTRANSLATE (yychar);
      yydestruct ("Cleanup: discarding lookahead",
                  yytoken, &yylval, &yylloc, context);
    }
  /* Do not reclaim the symbols of the rule which action triggered
     this YYABORT or YYACCEPT.  */
  YYPOPSTACK (yylen);
  YY_STACK_PRINT (yyss, yyssp);
  while (yyssp != yyss)
    {
      yydestruct ("Cleanup: popping",
		  yystos[*yyssp], yyvsp, yylsp, context);
      YYPOPSTACK (1);
    }
#ifndef yyoverflow
  if (yyss != yyssa)
    YYSTACK_FREE (yyss);
#endif
#if YYERROR_VERBOSE
  if (yymsg != yymsgbuf)
    YYSTACK_FREE (yymsg);
#endif
  /* Make sure YYID is used.  */
  return YYID (yyresult);
}



