#!/bin/sh

eval "$(/usr/bin/garden-exec)"
source $(dirname $0)/../env.sh
garden in-path PATH ${PREFIX:-`dirname $0`/../build}/bin
garden in-path PYTHONPATH ${PREFIX:-`dirname $0`/../build}/lib/python
garden load openeye-toolkits/2020.0.4-01c7/lib-python37
garden load rdkit/2019.03.1-05c7/lib-python37
garden load zstandard/0.11.1-01c7/lib-python37
garden load brotli/1.0.9-01c7/lib-python37
garden load valgrind/3.15.0-01c7/bin

garden with -m $PYTHON3/bin  pytest $(dirname $0) -p no:azurepipelines "$@"

