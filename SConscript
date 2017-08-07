Import('env')
import os
import sys
import subprocess

def _AddPython3Module( env, *args, **kwds ):
    ''' Install python files $args into lib/python3.
    Optional prefix keyword argument puts files into lib/python/$prefix/.
    '''
    prefix=kwds.get('prefix', '')
    progs=[]
    for s in args:
        py='$OBJDIR/lib/python3/%s/%s' % (prefix, os.path.basename(s))
        progs.extend( env.Command(py, s, [
            Copy("$TARGET", "$SOURCE"),
            'python3 -m compileall $TARGET']))
    _Install(env,_Join('lib/python3',prefix),progs)
    return progs


def _AddPythonExtension(env, name, *args, **kwds):
    env=env.Clone()
    incs=subprocess.check_output('python-config --includes'.split()).strip()
    incs += ' -I' + subprocess.check_output(['python', '-c', 'import numpy; print(numpy.get_include())']).strip()
    env.Prepend(CPPFLAGS=incs)
    kwds = kwds.copy()
    kwds.update(LIBPREFIX='', SHLIBSUFFIX='.so', SHOBJSUFFIX='.os')
    if env['PLATFORM']=='darwin':
      env.AppendUnique( LINKFLAGS=['-undefined','dynamic_lookup'] )
    prefix=kwds.get('prefix', '')
    py='$OBJDIR'+'/'+_Join('lib/python',prefix)+'/'+name
    lib=env.SharedLibrary(py, *args, **kwds)
    _Install(env,_Join('lib/python',prefix),lib)
    return lib


def _AddPython3Extension(env, name, *args, **kwds):
    env=env.Clone()
    incs=subprocess.check_output('python3-config --includes'.split()).strip()
    incs += ' -I' + subprocess.check_output(['python3', '-c', 'import numpy; print(numpy.get_include())']).strip()
    env.Prepend(CPPFLAGS=incs)
    kwds = kwds.copy()
    kwds.update(LIBPREFIX='', SHLIBSUFFIX='.so', SHOBJSUFFIX='.os3')
    if env['PLATFORM']=='darwin':
      env.AppendUnique( LINKFLAGS=['-undefined','dynamic_lookup'] )
    prefix=kwds.get('prefix', '')
    py='$OBJDIR'+'/'+_Join('lib/python3',prefix)+'/'+name
    lib=env.SharedLibrary(py, *args, **kwds)
    _Install(env,_Join('lib/python3',prefix),lib)
    return lib


def _Install(env,subdir,stuff,**kwds):
  if env['PREFIX'] is not None:
    tgt=_Join('$PREFIX',subdir)
    if kwds.get('SHLIBVERSION'):
        EnsureSConsVersion(2,3)
        out=env.InstallVersionedLib(tgt,stuff,**kwds)
    else:
        out=env.Install(tgt,stuff,**kwds)
    env.Alias('install', tgt)
    return out
  return None

def _Join(path,sub=''):
  if sub: return path + '/' + sub
  else:   return path

def Customize(env):
    for name, func in globals().items():
        if name.startswith('_Add') and callable(func):
            name=name[1:]
            env.AddMethod( func, name )


Customize(env)

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

env.AddShare('MODULES')

opts=Variables()
opts.Add("MSYS_WITHOUT_PYTHON", "without python?")
opts.Update(env)

if env.get("MSYS_WITHOUT_PYTHON"):
    print "MSYS_WITHOUT_PYTHON set; will not build python extensions or tools."
else:
    env.SConscript('python/SConscript')
    env.SConscript('tools/SConscript')

