#ifndef desres_msys_atomsel_expression_hxx
#define desres_msys_atomsel_expression_hxx

#include "predicate.hxx"
#include "keyword.hxx"

namespace desres { namespace msys { namespace atomsel {

  /* an expression holds a floating point value for each element in 
   * a selection. */
  struct Expression : public boost::enable_shared_from_this<Expression> {
    virtual ~Expression() {}
    virtual void eval( const Selection& s, std::vector<Dbl>& v ) = 0;
    virtual void dump(std::ostream& str) const = 0;
  };
  typedef boost::shared_ptr<Expression> ExpressionPtr;

  /* unary operand expression.  operand may be "+", "-", sqr, sqrt, abs. */
  ExpressionPtr unary_expression( const std::string& op, ExpressionPtr sub );

  /* binary operand expression.  Supports _, -, *, /, **, ^, and %. */
  ExpressionPtr binary_expression( const std::string& op, 
      ExpressionPtr lhs, ExpressionPtr rhs );

  /* an expression with the same literal value */
  ExpressionPtr literal_expression( double v );

  /* an expression whose value is a keyword.  Strings evaluate to 0;
   * booleans as 0 or 1. */
  ExpressionPtr keyword_expression( KeywordPtr key );

  enum RelationType {
    RELATION_EQ=0, 
    RELATION_NE, 
    RELATION_LT, 
    RELATION_LE, 
    RELATION_GE, 
    RELATION_GT
  };
  /* A predicate that compares two expressions */
  PredicatePtr relation_predicate( RelationType type,
      ExpressionPtr lhs, ExpressionPtr rhs );

}}} // ns

#endif
