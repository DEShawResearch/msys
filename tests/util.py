import os, sys
import platform

# DON'T CHANGE THIS! RUN ut.py out of the base directory.  -- JRG
desres_os  = os.getenv("DESRES_OS",  platform.system())
desres_isa = os.getenv("DESRES_ISA", platform.machine())
TMPDIR = os.getenv('TMPDIR') if platform.system() != 'Darwin' else None
TMPDIR = TMPDIR or 'objs/%s/%s' % (desres_os, desres_isa)
suffix = '3' if sys.version_info.major==3 else ''
sys.path.insert(0,os.path.join(TMPDIR, 'lib', 'python%s' % suffix))

