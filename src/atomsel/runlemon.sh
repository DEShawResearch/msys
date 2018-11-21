#!/bin/sh

cd `dirname $0`
lemon atomsel.y
mv atomsel.c atomsel.cxx

