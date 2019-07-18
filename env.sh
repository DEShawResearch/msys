PYTHON=desres-python/2.7.15-10c7
PYTHON3=desres-python/3.7.4-01c7
YAS=yas/0.143c7

garden prepend-path PYTHONPATH $(dirname $0)/../lib/python
garden load $PYTHON3/bin
garden load $YAS/lib-python37

