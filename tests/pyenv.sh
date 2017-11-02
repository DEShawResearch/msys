garden env-keep-only TMPDIR
source `dirname $0`/../MODULES
if [ "$1" == "-3" ]
then
   shift
   PYSUFFIX='3'
   garden load $PYTHON3/bin
else
   PYSUFFIX=''
   garden load $PYTHON/bin
fi
which python
