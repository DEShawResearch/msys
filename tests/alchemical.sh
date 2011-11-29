#!/usr/bin/env desres-exec
#{
# desres-with \
# -m ndiff/2.00-02/bin \
# -m maeff_inputs/1.4.1/share \
# -- sh $0 "$@"
#}

TMP=/tmp
set -e
DUMPOPTS="--without-provenance --reorder --without-forcefield"

for DIR in $DESMOND_MAEFF_INPUT_PATH/fep/*
#for DIR in $DESMOND_MAEFF_INPUT_PATH/fep/fep1
do
    echo "Reading A.mae, B.mae, and atom.map from $DIR"
    
    # viparr1 route.
    mae2dms $DIR/C.mae $TMP/viparr.dms
    sqlite3 /tmp/viparr.dms "update alchemical_dihedral_trig_term set p0=p3,p1=p2,p2=p1,p3=p0 where p1>p2"
    sqlite3 /tmp/viparr.dms "update alchemical_pair_12_6_es_term set p0=p1,p1=p0 where p0>p1"

    dms-dump $DUMPOPTS $TMP/viparr.dms | grep -v "\-----" > $TMP/viparr.dump
    
    # msys route.
    mae2dms $DIR/A.mae $TMP/A.dms
    mae2dms $DIR/B.mae $TMP/B.dms
    dms-alchemical $TMP/A.dms $TMP/B.dms $DIR/atom.map $TMP/msys.dms
    sqlite3 /tmp/msys.dms "update alchemical_dihedral_trig_term set p0=p3,p1=p2,p2=p1,p3=p0 where p1>p2"
    sqlite3 /tmp/msys.dms "update alchemical_pair_12_6_es_term set p0=p1,p1=p0 where p0>p1"

    # compare
    dms-dump $DUMPOPTS $TMP/msys.dms | grep -v "\-----" > $TMP/msys.dump
    ndiff -relerr 1e-4 -silent $TMP/viparr.dump $TMP/msys.dump

    echo "OK"
done

echo All tests passed.

