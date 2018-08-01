import os, sys
import platform

# DON'T CHANGE THIS! RUN ut.py out of the base directory.  -- JRG
TMPDIR = os.getenv('TMPDIR') if platform.system() != 'Darwin' else None
TMPDIR = TMPDIR or 'build'
sys.path.insert(0,os.path.join(TMPDIR, 'lib/python'))

