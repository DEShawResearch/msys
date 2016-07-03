%option reentrant
%option prefix="msys_smiles_"
%option bison-bridge
%option bison-locations
%option noyywrap
%option nounput

%{

#include "smiles/smiles.hxx"
#include "elements.hxx"
using namespace desres::msys::smiles;
using desres::msys::ElementForAbbreviation;

#include "parser.hxx"

#define YY_EXTRA_TYPE Smiles*

#define YY_INPUT(buf,result,max_size) \
{ \
    buf[0] = yyextra->getc(); \
    result = buf[0] ? 1 : YY_NULL; \
} \

%}

bond    [\-=#\$:/\\]
dot     [.]
percent [%]
digit   [0-9]
number  {digit}+

aliphatic_organic   (B|C|N|O|S|P|F|Cl|Br|I)
aromatic_organic    [bcnosp]
wild                "*"
atom                ({aliphatic_organic}|{aromatic_organic}|{wild})

element (H|He|Li|Be|B|C|N|O|F|Ne|Na|Mg|Al|Si|P|S|Cl|Ar|K|Ca|Sc|Ti|V|Cr|Mn|Fe|Co|Ni|Cu|Zn|Ga|Ge|As|Se|Br|Kr|Rb|Sr|Y|Zr|Nb|Mo|Tc|Ru|Rh|Pd|Ag|Cd|In|Sn|Sb|Te|I|Xe|Cs|Ba|Hf|Ta|W|Re|Os|Ir|Pt|Au|Hg|Tl|Pb|Bi|Po|At|Rn|Fr|Ra|Rf|Db|Sg|Bh|Hs|Mt|Ds|Rg|Cn|Fl|Lv|La|Ce|Pr|Nd|Pm|Sm|Eu|Gd|Tb|Dy|Ho|Er|Tm|Yb|Lu|Ac|Th|Pa|U|Np|Pu|Am|Cm|Bk|Cf|Es|Fm|Md|No|Lr)
aromatic ([bcnops]|se|as)
symbol  ({element}|{aromatic}|"*")

chiral  (@|@@)
hcount  H{digit}?
ncharge [\-+]{digit}?{digit}?
plus2   [+][+]
minus2  [\-][\-]
charge  ({plus2}|{minus2}|{ncharge})
class   :{number}

lbracket "["
rbracket "]"
%x bracket
%x symbol

%%

\(  { return LPAREN; }
\)  { return RPAREN; }

{atom}  { strcpy(yylval->sym, yytext); return ATOM; }

{lbracket}              { BEGIN(bracket); return LBRACKET; }

<bracket>{number}   { yylval->num = atoi(yytext); return ISOTOPE; }
<bracket>{symbol}   { strcpy(yylval->sym, yytext); BEGIN(symbol);return SYMBOL;}
<symbol>{chiral}    { yylval->num = strlen(yytext); return CHIRAL;}
<symbol>{hcount}    { yylval->num = yytext[1] ? atoi(yytext+1) : 1; return HCOUNT;}
<symbol>{class}     { yylval->num = atoi(yytext+1); return CLASS; }
<symbol>{charge}    { yylval->num = (yytext[1]=='-' ? -1 : 1) * (yytext[0]=='+' ? 1 : -1); return CHARGE; }

<symbol>{rbracket}  { BEGIN(INITIAL); return RBRACKET; }

{bond}{digit}  { 
                    yylval->ringbond = yyextra->makeRingbond(*yytext,atoi(yytext+1));
                    return RINGBOND;
               }

{bond}{percent}{digit}{digit}  { 
                    yylval->ringbond = yyextra->makeRingbond(*yytext,atoi(yytext+2));
                    return RINGBOND;
               }

{digit}         {   yylval->ringbond = yyextra->makeRingbond(0,      atoi(yytext));
                    return RINGBOND;
                }
{percent}{digit}{digit} {   
                    yylval->ringbond = yyextra->makeRingbond(0,      atoi(yytext+1));
                    return RINGBOND;
                }
{bond} { yylval->bond = *yytext; return SBOND; }
{dot}  { return DOT;                           }

<INITIAL,bracket,symbol>. { yyextra->error = "Bad character: "; 
                            yyextra->error+=yytext; 
                            return ERR; }

%%

void Smiles::init_scanner() {
    yylex_init(&scanner);
    yyset_extra(this, scanner);
}

void Smiles::destroy_scanner() {
    yylex_destroy(scanner);
}

