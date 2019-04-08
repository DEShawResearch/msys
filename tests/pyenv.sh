garden env-keep-only TMPDIR
source `dirname $0`/../env.sh
if [ "$1" == "-3" ]
then
   shift
   PYSUFFIX='3'
   garden load $PYTHON3/bin
   garden load $YAS/lib-python3
else
   PYSUFFIX=''
   garden load $PYTHON/bin
   garden load $YAS/lib-python
fi
which python
garden in-path PATH ${PREFIX:-`dirname $0`/../build}/bin

