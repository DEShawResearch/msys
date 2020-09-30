#!/usr/bin/env desres-exec
#{
# exec desres-cleanenv \
# -m msys/1.7.54/lib-python \
# -m dms_inputs/1.5.8/share \
# -- python $0 "$@"
#}

import sys, os

TMPDIR = os.getenv("TMPDIR", "objs/Linux/x86_64")
sys.path.insert(0, os.path.join(TMPDIR, "lib", "python"))
import pfx

import numpy as NP
import msys
from time import time


def bench(path, fixbonds=False, glue=[], center=None, fit=None):
    fullpath = "%s/%s" % (os.environ["DMS_INPUTS_PATH"].split(":")[0], path)
    mol = msys.Load(fullpath)
    box = mol.cell
    dpos = mol.positions
    fpos = dpos.astype("f")
    dpos2 = dpos + [1, 2, 3]
    fpos2 = dpos2.astype("f")
    for pos in fpos, dpos:
        print "%s data from %s" % ("single" if pos is fpos else "double", path)
        pos2 = fpos2 if pos is fpos else dpos2
        kls = pfx.Pfx
        t0 = time() * 1000
        p = kls(mol.topology, fixbonds=fixbonds)
        t1 = time() * 1000
        if center:
            p.align(mol.selectIds(center))
        elif fit:
            ids = mol.selectIds(fit)
            p.align(ids, pos2[ids])
        t2 = time() * 1000
        p.apply(pos, box)
        t3 = time() * 1000
        print "%s\tconstruct %8.3fms center %8.3fms apply %8.3fms" % (
            "single" if kls is pfx.Pfx else "double",
            t1 - t0,
            t2 - t1,
            t3 - t2,
        )


bench("apoa1.dms", fixbonds=True, glue="backbone", center="protein")
bench("small_vancomycin_complex.dms", glue="backbone")
bench("small_vancomycin_complex.dms", glue="backbone", center="residue 0")
bench("1vcc.dms", fit="backbone")
bench("1vcc.dms", center="backbone")
