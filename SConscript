Import('env')
import os
import sys
import subprocess

# Avoid doing dependency checks on garden files.
if True:
    cpp=[]
    flg=[]
    for p in env['CPPPATH']:
        if p.startswith('/proj') or p.startswith('/gdn'):
            flg.append('-I%s' % p)
        else:
            cpp.append(p)
    env.Replace(CPPPATH=cpp)
    env.Append(CFLAGS=flg, CXXFLAGS=flg)

if "SCHRODINGER_SRC" not in os.environ:
    env.Append(
            # SSE2 for src/within.hxx.  It's optional, but way way slower without.
            CCFLAGS='-O2 -g -msse4.1',
            CFLAGS='-Wall',
            # sadly, need -Wno-deprecated-declarations because of boost.
            CXXFLAGS='-std=c++14 -Wall -Werror -Wno-deprecated-declarations',
            CPPDEFINES=[
                'BOOST_SYSTEM_NO_DEPRECATED',
                ],
            LINKFLAGS='-g'
            )

if env['PLATFORM']=='darwin':
    env.Append(CXXFLAGS='-ftemplate-depth=500')

env.SConsignFile('%s/.sconsign' % (env['OBJDIR'].strip('#')))

env.SConscript('src/SConscript')
if "SCHRODINGER_SRC" not in os.environ:
    env.SConscript('tests/SConscript')

env.AddShare('env.sh')

opts=Variables()
opts.Add("MSYS_WITHOUT_PYTHON", "without python?")
opts.Update(env)

if env.get("MSYS_WITHOUT_PYTHON"):
    print "MSYS_WITHOUT_PYTHON set; will not build python extensions or tools."
else:
    env.SConscript('python/SConscript')
    env.SConscript('tools/SConscript')

