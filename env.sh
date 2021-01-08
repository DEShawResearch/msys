export PYTHON37=desres-python/3.7.7-06c7
PYTHON3=desres-python/3.9.1-01c7
YAS=yas/0.164c7

garden prepend-path PYTHONPATH $(dirname $0)/../lib/python
garden load $PYTHON3/bin
garden load $YAS/lib-python37

