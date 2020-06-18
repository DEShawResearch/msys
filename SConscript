Import('env')
import os
import sys
import subprocess

#
# Inherit environment from outside for these environment variables
#
TOOLCHAIN = ["CC", "CXX", "LD", "AR", "AS", "STRIP"]
for envvar in TOOLCHAIN:
    if envvar in os.environ:
        env.Replace(**{envvar: os.environ[envvar]})

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

env.Append(
        # SSE2 for src/within.hxx.  It's optional, but way way slower without.
        CCFLAGS=['-O2', '-g', '-msse4.1'],
        CFLAGS='-Wall',
        # sadly, need -Wno-deprecated-declarations because of boost.
        CXXFLAGS="-std=c++11 -Wall -Wno-deprecated-declarations",
        CPPDEFINES=[
            'BOOST_SYSTEM_NO_DEPRECATED',
            ],
        LINKFLAGS='-g'
        )

if env['PLATFORM']=='darwin':
    env.Append(CXXFLAGS='-ftemplate-depth=500')

env.SConsignFile('%s/.sconsign' % (env['OBJDIR'].strip('#')))

env.SConscript('src/SConscript')
env.SConscript('tests/SConscript')

env.SConscript('external/inchi/SConscript')
env.SConscript('external/lpsolve/SConscript')

env.AddShare('env.sh')
env.SConscript('python/SConscript')
env.SConscript('tools/SConscript')

if 'BUILD_WHEEL' in os.environ:
    env['WHEEL_DIR'] = 'wheel'
    env.AddWheel('pyproject.toml', pyver=os.getenv('BUILD_WHEEL_VERSION'))

