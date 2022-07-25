PYTHON3=desres-python/3.10.5dev1-04c7
YAS=yas/0.189c7
# current yas, as of, is yas/0.189c7

garden prepend-path PYTHONPATH $(dirname $0)/../lib/python
garden load $PYTHON3/bin
garden load $YAS/lib-python37
garden load $YAS/lib-python310
