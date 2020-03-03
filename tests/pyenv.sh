garden env-keep-only TMPDIR PREFIX
source `dirname $0`/../env.sh
garden load $PYTHON3/bin
garden load $YAS/lib-python37
garden in-path PATH ${PREFIX:-`dirname $0`/../build}/bin

