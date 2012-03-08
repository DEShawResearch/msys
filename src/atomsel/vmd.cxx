#include "vmd.hxx"
#include "msys_keyword.hxx"
#include "vmd_keyword.hxx"
#include "keyword_predicate.hxx"

#include "keyword.hxx"
#include "expression.hxx"
#include "within_predicate.hxx"
#include "VmdLexer.h"
#include "VmdParser.h"

#include <stdexcept>
#include <sstream>

using namespace desres::msys::atomsel;
using namespace desres::msys::atomsel::vmd;
using desres::msys::SystemPtr;

#define THROW_FAILURE( args ) do { \
    std::stringstream ss; \
    ss << args; \
    throw std::runtime_error(ss.str()); \
} while (0)


namespace {
  typedef ANTLR3_BASE_TREE Tree;
  PredicatePtr  parse( Tree * tree, SystemPtr sys, StrList& s );

  static bool contains(StrList& v, const char *s) {
      for (unsigned i=0; i<v.size(); i++) {
          if (!strcmp(v[i], s)) return true;
      }
      return false;
  }
}

PredicatePtr desres::msys::atomsel::vmd::parse(
        const std::string& sel, SystemPtr sys, StrList& prev) {

  char *txt = const_cast<char*>(sel.c_str());

  pANTLR3_INPUT_STREAM input=antlr3NewAsciiStringInPlaceStream(
      (pANTLR3_UINT8)txt, sel.size(), (pANTLR3_UINT8)"selection");

  pVmdLexer lxr = desres_msys_VmdLexerNew(input);

  pANTLR3_COMMON_TOKEN_STREAM tstream = antlr3CommonTokenStreamSourceNew(
      ANTLR3_SIZE_HINT, TOKENSOURCE(lxr));

  pVmdParser psr = desres_msys_VmdParserNew(tstream);

  VmdParser_start_return ast=psr->start(psr);

  if (psr->pParser->rec->state->errorCount > 0) {
    int pos=psr->pParser->rec->state->exception->charPositionInLine;
    std::stringstream ss;
    ss << "VmdGrammar: parsing failed at offset " 
        << pos 
        << "around '" << sel.substr(pos, std::string::npos).c_str() << "'"
        ;
    throw std::runtime_error(ss.str());
  }

  PredicatePtr pred = ::parse(ast.tree, sys, prev);

  /* TODO: create scoped_ptr for parser state in case exception is thrown
   * during parsing. */
  psr->free(psr);
  tstream->free(tstream);
  lxr->free(lxr);
  input->close(input);

  return pred;
}

namespace {


  const char* find_macro(std::string const& name, SystemPtr mol) { 
      std::string const& m = mol->selectionMacroDefinition(name);
      return m.size() ? m.c_str() : NULL;
  }

  char *str(Tree *tree) {
    static char empty[]="";
    if (!tree) return empty;
    pANTLR3_STRING s=tree->getText(tree);
    return s ? reinterpret_cast<char*>(s->chars) : empty;
  }

  Tree *child(Tree *tree, unsigned int i) {
    if (!tree) return NULL;
    if (!tree->children) return NULL;
    if (tree->children->count < i) return NULL;
    return (Tree*)tree->children->get(tree->children, i);
  }
  PredicatePtr parseAnd( Tree *tree, SystemPtr ent, StrList& s ) {
    return and_predicate( 
        parse(child(tree,0), ent, s),
        parse(child(tree,1), ent, s) );
  }

  PredicatePtr parseOr( Tree *tree, SystemPtr ent, StrList& s ) {
    return or_predicate(
        parse(child(tree,0), ent, s),
        parse(child(tree,1), ent, s) );
  }

  PredicatePtr parseNot( Tree *tree, SystemPtr ent, StrList& s ) {
    return not_predicate(
        parse(child(tree,0), ent, s) );
  }

  KeywordPtr get_keyword( const char * id, SystemPtr ent ) {
    KeywordPtr key = 
      !strcmp(id,"atomicnumber") ? keyword_anum(ent) :
      !strcmp(id,"chain") ? keyword_chain(ent) :
      !strcmp(id,"charge") ? keyword_charge(ent) :
      !strcmp(id,"fragment") ? keyword_fragment(ent) :
      !strcmp(id,"index") ? keyword_index(ent) :
      !strcmp(id,"mass") ? keyword_mass(ent) :
      !strcmp(id,"name") ? keyword_name(ent) :
      !strcmp(id,"numbonds") ? keyword_numbonds(ent) :
      !strcmp(id,"resid") ? keyword_resid(ent) :
      !strcmp(id,"residue") ? keyword_residue(ent) :
      !strcmp(id,"resname") ? keyword_resname(ent) :
      !strcmp(id,"fragid") ? keyword_fragid(ent) :
      !strcmp(id,"x") ? keyword_x(ent) :
      !strcmp(id,"y") ? keyword_y(ent) :
      !strcmp(id,"z") ? keyword_z(ent) :
      !strcmp(id,"vx") ? keyword_x(ent) :
      !strcmp(id,"vy") ? keyword_y(ent) :
      !strcmp(id,"vz") ? keyword_z(ent) :
      KeywordPtr();

    if (!key) key = keyword_atomprop(ent, id);
    return key;
  }

  KeywordPtr get_boolean( const char * id, SystemPtr ent ) {
    KeywordPtr key = 
      !strcmp(id,"water") ? keyword_water(ent) :
      !strcmp(id,"hydrogen") ? keyword_hydrogen(ent) :
      !strcmp(id,"backbone") ? keyword_backbone(ent) :
      !strcmp(id,"protein") ? keyword_protein(ent) :
      //!strcmp(id,"alchemical") ? keyword_alchemical(ent) :
      !strcmp(id,"nucleic") ? keyword_nucleic(ent) :
      KeywordPtr();
    return key;
  }

  PredicatePtr parseSame( Tree *tree, SystemPtr ent, StrList& s ) {
    const char *id = str(child(tree,0));
    KeywordPtr key = get_keyword(id,ent);
    if (!key) {
      THROW_FAILURE(
          "VmdGrammar: unrecognized keyword '" << std::string(id)
          << "'");
    }
    return same_keyword_predicate( key, full_selection(ent), 
            parse(child(tree,1), ent, s) );
  }

  PredicatePtr parseKeyword( Tree *tree, SystemPtr ent, StrList& s ) {
    const char *id = str(child(tree,0));
    int nterms = tree->children->count-1;
    KeywordPtr key = get_keyword(id,ent);
    if (key) {
      if (nterms==0) {
        THROW_FAILURE("keyword '" << id << "' was not followed by any terms");
      }
      boost::shared_ptr<KeywordPredicate> pred(new KeywordPredicate(key));
      for (int i=0; i<nterms; i++) {
        Tree * c = child(tree,i+1);
        int type=c->getType(c);
        if (type==TO) { /* range */
          pred->addRange(Range(str(child(c,0)),str(child(c,1))));

        } else if (type==LIT) {
          pred->addLiteral(str(c));

        } else if (type==REGEX) { 
          std::string r(str(child(c,0))+1);
          r.resize(r.size()-1);
          pred->addRegex(Regex(r));

        } else {
          THROW_FAILURE("Unexpected type " << type << " following key " << id);
        }
      }
      return pred;
    }
    key = get_boolean(id, ent);
    if (key) {
      if (nterms!=0) {
        THROW_FAILURE("boolean '" << id << "' followed by terms");
      }
      return boolean_predicate(key);
    }

    /* special case "all" and "none" */
    if (!strcmp(id,"all")) return all_predicate();
    if (!strcmp(id,"none")) return none_predicate();

    /* check for macros defined in the system */
    const char* user_macro = find_macro(id, ent);
    if (user_macro) {
        if (contains(s,user_macro)) {
            MSYS_FAIL("Recursion detected in selection macro " << id);
        }
        s.push_back(user_macro);
        PredicatePtr p = desres::msys::atomsel::vmd::parse(user_macro,ent,s);
        s.pop_back();
        return p;
    }

    THROW_FAILURE("VmdGrammar: unrecognized keyword " << id);
    return PredicatePtr();
  }

  /* forward declaration */
  ExpressionPtr parseExpression(Tree* tree, SystemPtr ent);

  ExpressionPtr parseUnaryExpression(Tree* tree, SystemPtr ent) {
    const char *id = str(child(tree,0));
    ExpressionPtr sub = parseExpression(child(tree,1), ent);
    return unary_expression(id,sub);
  }

  ExpressionPtr parseBinaryExpression(Tree* tree, SystemPtr ent) {
    const char *id = str(child(tree,0));
    ExpressionPtr sub1 = parseExpression(child(tree,1), ent);
    ExpressionPtr sub2 = parseExpression(child(tree,2), ent);
    return binary_expression(id,sub1,sub2);
  }

  ExpressionPtr parseLiteralExpression(Tree * tree, SystemPtr ent ) {
    const char * tok = str(child(tree,0));
    KeywordPtr key = get_keyword(tok,ent);
    ExpressionPtr expr;
    if (key) { 
      expr = keyword_expression(key);
    } else {
      char * endptr;
      double val = strtod(tok, &endptr);
      if (endptr!=tok) {
        expr = literal_expression(val);
      } else {
        THROW_FAILURE("Could not understand expression '" << tok << "'");
      }
    }
    return expr;
  }

  ExpressionPtr parseExpression(Tree* tree, SystemPtr ent) {
    switch (tree->getType(tree)) {
      case BINARYOP: return parseBinaryExpression(tree,ent); break;
      case LITERAL: return parseLiteralExpression(tree,ent); break;
      case FUNCTION: return parseUnaryExpression(tree,ent); break;
      case UNARYOP: return parseUnaryExpression(tree,ent); break;
    }
    THROW_FAILURE("Unsupported expression type " << tree->getType(tree));
  }

  PredicatePtr parseRelation(Tree* tree, SystemPtr ent) {
    Tree *op = child(tree,0);
    int type=op->getType(op);

    RelationType c;
    switch (type) {
      case LESS:   c=RELATION_LT; break;
      case MORE:   c=RELATION_GT; break;
      case LESSEQ: c=RELATION_LE; break;
      case MOREEQ: c=RELATION_GE; break;
      case NEQUAL: c=RELATION_NE; break;
      default:     c=RELATION_EQ;
    }

    return relation_predicate(c, parseExpression(child(tree,1),ent),
                                 parseExpression(child(tree,2),ent) );
  }

  PredicatePtr parseWithin(Tree* tree, SystemPtr ent, StrList& s) {
    double rad = atof(str(child(tree,0)));
    return within_predicate(ent,rad,parse(child(tree,1),ent,s));
  }

  PredicatePtr parseExwithin(Tree* tree, SystemPtr ent, StrList& s) {
    double rad = atof(str(child(tree,0)));
    return exwithin_predicate(ent,rad,parse(child(tree,1),ent,s));
  }

  PredicatePtr parsePbwithin(Tree* tree, SystemPtr ent, StrList& s) {
    double rad = atof(str(child(tree,0)));
    return pbwithin_predicate(ent,rad,parse(child(tree,1),ent,s));
  }

  PredicatePtr parseWithinBonds(Tree* tree, SystemPtr ent, StrList& s) {
    int n = atoi(str(child(tree,0)));
    return withinbonds_predicate(ent,n,parse(child(tree,1),ent,s));
  }

  PredicatePtr parse( Tree *tree, SystemPtr ent, StrList& s ) {
    if (!tree) THROW_FAILURE("atomsel::vmd - Unexpected NULL tree");
    switch (tree->getType(tree)) {
      case AND: return parseAnd(tree,ent,s);
      case OR:  return parseOr(tree,ent,s);
      case NOT: return parseNot(tree,ent,s);
      case KEYWORD: return parseKeyword(tree,ent,s);
      case SAME: return parseSame(tree,ent,s);
      case RELATION: return parseRelation(tree,ent);
      case WITHIN: return parseWithin(tree,ent,s);
      case EXWITHIN: return parseExwithin(tree,ent,s);
      case PBWITHIN: return parsePbwithin(tree,ent,s);
      case WITHINBONDS: return parseWithinBonds(tree,ent,s);
      default:;
    }
    THROW_FAILURE(
            "VmdGrammar parser: unexpected node type " << tree->getType(tree));
  }
}

