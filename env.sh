PYTHON3=desres-python/3.7.7-06c7
YAS=yas/0.161c7

garden prepend-path PYTHONPATH $(dirname $0)/../lib/python
garden load $PYTHON3/bin
garden load $YAS/lib-python37

