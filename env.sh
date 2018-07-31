PYTHON=desres-python/2.7.15-08c7
PYTHON3=desres-python/3.6.6-02c7
YAS=yas/0.111-beta-c7

garden prepend-path PYTHONPATH $(dirname $0)/../python
garden load $PYTHON/bin
garden load $YAS/lib-python

