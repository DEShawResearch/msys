#!/usr/bin/env desres-exec
#{
# exec desres-cleanenv \
# -m flex/2.5.35-1/bin \
# -m bison/2.5-11A/bin \
# -- sh $0 "$@"
#}

flex -o lexer.cxx vmd.l
bison -o parser.cxx vmd.y
