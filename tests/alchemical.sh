#!/usr/bin/env desres-exec
#{
# exec desres-cleanenv \
# -m ndiff/2.00-02/bin \
# -m dms_inputs/1.5.8/share \
# -m msys/1.7.84/bin \
# -- sh $0 "$@"
#}

TMP=/tmp
set -e
DUMPOPTS="--without-provenance --reorder --without-forcefield --without-paraminfo"

for DIR in $DMS_INPUTS_PATH/fep/fep*
#for DIR in $DESMOND_MAEFF_INPUT_PATH/fep/fep1
do
    rm -f old.dms new.dms

    echo "Reading A.dms, B.dms, and atom.map from $DIR"
    
    # old msys route.
    dms-alchemical $DIR/A.dms $DIR/B.dms $DIR/atom.map $TMP/old.dms
    dms-dump $DUMPOPTS $TMP/old.dms | grep -v "\-----" > $TMP/old.dump
    
    # new msys route.
    ./objs/$DESRES_OS/$DESRES_ISA/bin/dms-alchemical $DIR/A.dms $DIR/B.dms $DIR/atom.map $TMP/new.dms
    dms-dump $DUMPOPTS $TMP/new.dms | grep -v "\-----" > $TMP/new.dump

    # compare
    ndiff -relerr 1e-4 -silent $TMP/old.dump $TMP/new.dump

    echo "OK"
done

echo All tests passed.

