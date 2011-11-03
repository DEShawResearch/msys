/* @COPYRIGHT@ */


grammar Vmd;

options {
  language = C;
  output   = AST;
  // backtracking is necessary in order to distinguish between a grouping 
  // of selections, and a group of expressions; e.g.,
  //   (backbone or resname GLY) and resid 10 to 20
  // vs
  //   (x+y)*z < 10
  backtrack = true;
}

tokens {
  AND = 'and' ;
  OR  = 'or'  ;
  LESS   = '<' ;
  LESSEQ = '<=';
  MORE   = '>' ;
  MOREEQ = '>=';
  EQUAL  = '==';
  NEQUAL = '!=';
  KEYWORD;
  FUNCTION;
  RELATION;
  UNARYOP;
  BINARYOP;
  LITERAL;
  REGEX;
  WITHIN = 'within';
  WITHINBONDS = 'withinbonds';
  EXWITHIN = 'exwithin';
  PBWITHIN = 'pbwithin';
  OF= 'of';
  SAME = 'same';
  AS = 'as';
  NOT = 'not';
  TO = 'to';
}

// The symbol table used to recognize keyword selections, booleans,
// and functions.
scope Symbols 
{
}

start 
  : orSelection EOF
    -> orSelection
  ;

orSelection 
  // a or b
  : (a=andSelection->$a) (OR b=andSelection-> ^(OR $orSelection $b) )*
  ;

andSelection
  // a and b
  : (a=selection->$a) (AND b=selection-> ^(AND $andSelection $b) )*
  ;

selection
  // ( selection )
  : '(' orSelection ')' 
    -> orSelection

  // not hydrogen
  | NOT selection
    -> ^(NOT selection)

  // same residue as selection 
  | SAME LIT AS orSelection
    -> ^(SAME LIT orSelection)

  // within 1.5 of protein
  | WITHIN LIT OF orSelection
    -> ^(WITHIN LIT orSelection)

  // exwithin 1.5 of protein
  | EXWITHIN LIT OF orSelection
    -> ^(EXWITHIN LIT orSelection)

  // pbwithin 1.5 of protein
  | PBWITHIN LIT OF orSelection
    -> ^(PBWITHIN LIT orSelection)

  // withinbonds 3 of residue 10 
  | WITHINBONDS LIT OF orSelection
    -> ^(WITHINBONDS LIT orSelection)

  // x < 5
  | (a=addExpr->$a) relOp (b=addExpr->$b) 
    -> ^(RELATION relOp $a $b)

  // resid 10 20 30 to 40
  | LIT literals
    -> ^(KEYWORD LIT literals) 
  ;

literals
  : ( LIT | range | regex ) *
  ;

range
  // 3 to 5
  : LIT TO LIT
    -> ^(TO LIT LIT)
  ;

regex  
  : REGEXVAL 
    -> ^(REGEX REGEXVAL)
  ;

relOp
  : LESS | MORE | LESSEQ | MOREEQ | EQUAL | NEQUAL
  ;

addExpr
  : (a=multExpr->$a) ( addOp b=multExpr-> ^(BINARYOP addOp $addExpr $b) )* 
  ;

addOp
  : '+' | '-'
  ;

multExpr
  : (a=unaryExpr->$a) ( multOp b=unaryExpr-> ^(BINARYOP multOp $multExpr $b) )* 
  ;

multOp
  : '*' | '/' 
  ;

unaryExpr
  : unaryOp expExpr -> ^(UNARYOP unaryOp expExpr )
  | expExpr
  ;

unaryOp
  : '+' | '-'
  ;

expExpr
  : (a=factor->$a) ( expOp b=factor-> ^(BINARYOP expOp $expExpr $b) )* 
  ;

expOp
  : '**' | '^'
  ;

factor
  : '(' addExpr ')' 
    -> addExpr

  | LIT '(' addExpr ')' 
    -> ^(FUNCTION LIT addExpr)

  | LIT 
    -> ^(LITERAL LIT)
  ;

LIT
  // we don't try hard at all to catch errors in literal values.
  : ( ALPHA | DIGIT | '.' | '_' | '\'' | '*' ) +
  ;
  
REGEXVAL
  //: (LIT | '[' | ']' |'+' | '-' )+
  //:	(~('\"'|' '|'\n' | '\r'))+
   :    '\"' (~'\"')* '\"'
  ;

fragment ALPHA : 'a'..'z' | 'A'..'Z';
fragment DIGIT : '0'..'9' ;
WS          : (' '|'\n' | '\r')    {$channel=HIDDEN;} ;

