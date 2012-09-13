#!/usr/bin/env desres-exec
#{
# exec desres-cleanenv \
# -m ndiff/2.00-02/bin \
# -m dms_inputs/1.5.8/share \
# -m msys/1.6.2/bin \
# -- sh $0 "$@"
#}

TMP=/tmp
set -e
DUMPOPTS="--without-provenance --reorder --without-forcefield --without-paraminfo"

canonicalize() {
    set +e
    sqlite3 $1 "update stretch_harm_term set p0=p1,p1=p0 where p0>p1"
    sqlite3 $1 "update angle_harm_term set p0=p2,p2=p0 where p0>p2"
    sqlite3 $1 "update dihedral_trig_term set p0=p3,p1=p2,p2=p1,p3=p0 where p1>p2"
    sqlite3 $1 "update pair_12_6_es_term set p0=p1,p1=p0 where p0>p1"
    sqlite3 $1 "update alchemical_stretch_harm_term set p0=p1,p1=p0 where p0>p1"
    sqlite3 $1 "update alchemical_angle_harm_term set p0=p2,p2=p0 where p0>p2"
    sqlite3 $1 "update alchemical_dihedral_trig_term set p0=p3,p1=p2,p2=p1,p3=p0 where p1>p2"
    sqlite3 $1 "update alchemical_pair_12_6_es_term set p0=p1,p1=p0 where p0>p1"
    set -e
}

for DIR in $DMS_INPUTS_PATH/fep/fep*
#for DIR in $DESMOND_MAEFF_INPUT_PATH/fep/fep1
do
    rm -f old.dms new.dms

    echo "Reading A.dms, B.dms, and atom.map from $DIR"
    
    # old msys route.
    dms-alchemical $DIR/A.dms $DIR/B.dms $DIR/atom.map $TMP/old.dms
    canonicalize $TMP/old.dms
    dms-dump $DUMPOPTS $TMP/old.dms | grep -v "\-----" > $TMP/old.dump
    
    # new msys route.
    ./objs/Linux/x86_64/bin/dms-alchemical --keep-alchemical-noop $DIR/A.dms $DIR/B.dms $DIR/atom.map $TMP/new.dms
    canonicalize $TMP/new.dms
    dms-dump $DUMPOPTS $TMP/new.dms | grep -v "\-----" > $TMP/new.dump

    # compare
    ndiff -relerr 1e-4 -silent $TMP/old.dump $TMP/new.dump

    echo "OK"
done

echo All tests passed.

