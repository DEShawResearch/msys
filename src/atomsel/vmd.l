%option reentrant
%option prefix="vmd_"
%option bison-bridge
%option bison-locations
%option noyywrap
%option nounput

%{

#include "vmd.hxx"
using namespace desres::msys::atomsel;

#include "parser.hxx"


#define YY_EXTRA_TYPE   VMD*

#define YY_INPUT(buf,result,max_size) \
{ \
    buf[0] = yyextra->getc(); \
    result = buf[0] ? 1 : YY_NULL; \
} \

%}

float           ([0-9]+\.|([0-9]*\.[0-9]+)([eE][-+]?[0-9]+)?)
digit           [0-9]
alpha           [a-zA-Z_]
alphanum        ({alpha}|{digit})
integer         ([-+]?){digit}+

weird           [\'\"]
like_1A         {digit}+({alpha}|{weird})+({alphanum}|{weird})*
like_C5prime    {alpha}+{alphanum}*{weird}+({alphanum}|{weird})*

%%

\(                      { return '('; }
\)                      { return ')'; }
(and|\&\&)              { return AND; }
(or|\&\&)               { return OR;  }
not                     { return NOT; }
within                  { return WITHIN; }
exwithin                { return EXWITHIN; }
pbwithin                { return PBWITHIN; }
withinbonds             { return WITHINBONDS; }
nearest                 { return NEAREST; }
pbnearest               { return PBNEAREST; }
same                    { return SAME; }
as                      { return AS; }
of                      { return OF;  }
(to|\.\.\.)             { return TO;  }


\"([^\"]|\\\\|\\\")*\"  { /* escapes \ and " in a "string"  */
                            /* remove initial and trailing quotes */
                            char* s = (char *)malloc(yyleng);
                            char* pat = s;
                            for (char* p=yytext+1; *p; ++p, ++s) {
                              /* replace escaped quote with quote */
                              if (p[0]=='\\' && (p[1]=='"' || p[1]=='\\')) ++p;
                              *s = *p;
                            }
                            s[-1]=0;
                            unsigned flags = Regex::perl;
                            try {
                                yylval->regex = new Regex(pat,flags);
                            } catch (std::exception& e) {
                                yyextra->error = e.what();
                                return ERR;
                            }
                            return REGEX;
                        }

\'([^\']|\\\\|\\\')*\'  { /* escapes \ and ' in a 'string', */
                            yylval->sval = (char *)malloc(yyleng);
                            char* s = yylval->sval;
                            for (char* p=yytext+1; *p; ++p, ++s) {
                              /* replace escaped quote with quote */
                              if (p[0]=='\\' && (p[1]=='\'' || p[1]=='\\')) ++p;
                              *s = *p;
                            }
                            s[-1]=0;
                            return SVAL;
                        }


{integer}           { try {
                        yylval->ival = parse_int(yytext); 
                        return IVAL;
                      } catch (std::exception& e) {
                        yyextra->error = e.what();
                        return ERR;
                      }
                    }


\"[^\"]*    { return ERR; } /* unterminated double-quoted string */
\'[^\']*    { return ERR; } /* unterminated single-quoted string */

"=="                { yylval->ival = RELATION_EQ; return CMP; }
"!="                { yylval->ival = RELATION_NE; return CMP; }
"<"                 { yylval->ival = RELATION_LT; return CMP; }
"<="                { yylval->ival = RELATION_LE; return CMP; }
">="                { yylval->ival = RELATION_GE; return CMP; }
">"                 { yylval->ival = RELATION_GT; return CMP; }
"+"                 { return ADD; }
"-"                 { return SUB; }
"*"                 { return MUL; }
"/"                 { return DIV; }
"%"                 { return MOD; }
"^"                 { return EXP; }
"**"                { return EXP; }
                          
[ \t\n]             ;

{like_1A}|{like_C5prime} { /* catch raw names like: 5' 1A C4' 3'A */
                           yylval->sval = strdup(yytext);
                           return SVAL;
                        }

{float}             { yylval->fval = atof(yytext); return FVAL; }


({alpha}{alphanum}*)    { 
                          yylval->key = yyextra->find_single(yytext);
                          if (yylval->key) return SINGLE;
                          yylval->key = yyextra->find_key(yytext);
                          if (yylval->key) 
                              return yylval->key->type == KEY_INT ? INTKEY :
                                     yylval->key->type == KEY_DBL ? FLTKEY :
                                                                    STRKEY;
                          if (yyextra->find_function(yytext)) {
                            yylval->sval = strdup(yytext);
                            return FUNC;
                          }
                          yylval->pred = yyextra->find_macro(yytext);
                          if (yylval->pred) return MACRO;

                          yylval->pred = yyextra->find_strfctn(yytext);
                          if (yylval->pred) return STRFCTN;
                          yylval->sval = strdup(yytext); 
                          return SVAL; 
                        }

.                       { yyextra->error = "Bad character: "; yyextra->error += yytext;
                          return ERR;
                        }
%%

void VMD::init_scanner() {
    yylex_init(&scanner);
    yyset_extra(this, scanner);
}

void VMD::destroy_scanner() {
    yylex_destroy(scanner);
}

