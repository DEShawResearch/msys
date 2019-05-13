PYTHON=desres-python/2.7.15-10c7
PYTHON3=desres-python/3.7.3-03c7
YAS=yas/0.139c7

garden prepend-path PYTHONPATH $(dirname $0)/../lib/python
garden load $PYTHON3/bin
garden load $YAS/lib-python37

