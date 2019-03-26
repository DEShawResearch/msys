PYTHON=desres-python/2.7.15-08c7
PYTHON3=desres-python/3.7.2-08c7
YAS=yas/0.130c7

garden prepend-path PYTHONPATH $(dirname $0)/../lib/python
garden load $PYTHON/bin
garden load $YAS/lib-python

