#!/usr/bin/env python

version = "1.7.323"
major_version, minor_version, micro_version = [int(_x) for _x in version.split('.')]
hexversion = (major_version << 16) | \
             (minor_version << 8)  | \
             (micro_version << 0)

if __name__ == '__main__':
    print(version)
