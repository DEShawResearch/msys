#!/bin/sh

cd $(dirname $0)
flex -o lexer.cxx smiles.l
bison -t --locations -o parser.cxx smiles.y

