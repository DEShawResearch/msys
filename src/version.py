#!/usr/bin/env python

from __future__ import print_function

version = "1.7.262"
major_version, minor_version, micro_version = [int(_x) for _x in version.split('.')]
hexversion = (major_version << 16) | \
             (minor_version << 8)  | \
             (micro_version << 0)

if __name__ == '__main__':
    print(version)

