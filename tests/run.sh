#!/bin/sh

eval "$(/usr/bin/garden-exec)"
source $(dirname $0)/../env.sh
garden in-path PATH ${PREFIX:-`dirname $0`/../build}/bin
garden in-path PYTHONPATH ${PREFIX:-`dirname $0`/../build}/lib/python
garden load openeye-toolkits/2019.4.2-02c7/lib-python37
garden load rdkit/2019.03.1-05c7/lib-python37
exec pytest $(dirname $0) "$@"

