#!/usr/bin/env python2.7

import sys
import os

TMPDIR = os.getenv("TMPDIR", "objs/Linux/x86_64")
sys.path.insert(0, os.path.join(TMPDIR, "lib", "python"))
import msys
import json


def selections():
    for line in open("tests/atomsel_inputs.txt"):
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        print line
        yield line


cfg = "tests/atomsel_tests.json"
with open(cfg) as f:
    d = json.loads(f.read())
paths = map(str, d.keys())

for p in paths:
    mol = msys.Load(p)
    results = [[s, mol.selectIds(s)] for s in selections()]
    d[p] = results

with open(cfg, "w") as f:
    json.dump(d, f, separators=(", \n", ": "))
