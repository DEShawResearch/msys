'''
DESRES extensions for SCons
'''
from __future__ import print_function

from SCons.Script import *
import os
import subprocess

EnsurePythonVersion(2,7)
EnsureSConsVersion(2,4)

_wheel_targets = dict(platlib=[], purelib=[], data=[], scripts=[])

def _AddObject(env, *args, **kwds):
    objdir = env.subst('$OBJDIR')
    root = Dir('#').path
    suffix = kwds.get('SHOBJSUFFIX') or env['SHOBJSUFFIX']
    objs = list()
    for node in env.arg2nodes(args):
        if node.get_suffix().startswith('.o'):
            objs.append(node)
            continue
        rel = os.path.relpath(node.srcnode().get_path(), root)
        name = os.path.splitext(rel)[0]
        dst = '%s/%s%s' % (objdir, name, suffix)
        obj = env.SharedObject(dst, node, **kwds)
        objs.append(obj)
    return objs

def _AddLibrary(env, name, *args, **kwds):
    objs = env.AddObject(*args, **kwds)
    if not isinstance(name, str):
        raise ValueError("AddLibrary: expected name to be string, got %s" % type(name))
    if 'archive' in kwds and kwds['archive']:
        lib=env.Library('$LIBDIR/%s' % name, objs, **kwds)
    else:
        if env['PLATFORM']=='darwin' and env.get('PREFIX'):
            env=env.Clone()
            env.Append(LINKFLAGS=['-install_name', '$PREFIX/lib/lib%s$SHLIBSUFFIX' % name])
        lib=env.SharedLibrary('$LIBDIR/%s' % name, objs, **kwds)
    _Install(env,'lib',lib, **kwds)
    return lib

def _AddProgram(env, name, *args, **kwds):
    objs = env.AddObject(*args, **kwds)
    prog=env.Program('$BINDIR/%s' % name, objs, **kwds)
    _Install(env, 'bin', prog)
    wheel_prefix = 'bin'
    _wheel_targets['scripts'].append((prog, wheel_prefix))
    return prog

def add_test(env, prog):
    alias = env.Alias('test', [prog], prog[0].path)
    AlwaysBuild(alias)
    return prog

def _AddTestProgram(env, name, *args, **kwds):
    objs = env.AddObject(*args, **kwds)
    prog=env.Program('$TESTDIR/%s' % name, objs, **kwds)
    return add_test(env, prog)

def munge_header(env, source, target):
    dst = target[0].get_path()
    with open(dst) as fp:
        lines = fp.readlines()
    with open(dst, 'w') as fp:
        print('#!/usr/bin/env python', file=fp)
        print('from __future__ import print_function', file=fp)
        print('import os, sys', file=fp)
        print('sys.path.insert(0,os.path.dirname(__file__)+"/../lib/python")', file=fp)

        skip = lines[1].startswith('#{')
        for line in lines[1:]:
            if skip:
                if line.startswith('#}'):
                    skip = False
            else:
                fp.write(line)

def _AddScript( env, name, src=None ):
    if src is None: src=name
    name=os.path.basename(name)
    cmdlist = [
        Delete("$TARGET"),
        Copy("$TARGET", "$SOURCE"),
        Chmod("$TARGET", 0o755)]
    if not env['ENV'].get('DESRES_LOCATION'):
        cmdlist.append(munge_header)
    prog = env.Command('$BINDIR/%s'%name, src, cmdlist)
    _InstallAs(env,'bin'+'/'+name,prog)
    wheel_prefix = 'bin'
    _wheel_targets['scripts'].append((prog, wheel_prefix))
    return prog 

def _AddTestScript( env, name, src=None ):
    if src is None: src=name
    name=os.path.basename(str(name))
    prog=env.Command('$TESTDIR/%s'%name, src, [
        Delete("$TARGET"),
        Copy("$TARGET", "$SOURCE"),
        Chmod("$TARGET", 0o755)])
    return add_test(env, prog)

def _AddExampleProgram(env, name, *args, **kwds):
    objs = env.AddObject(*args, **kwds)
    prog=env.Program('$TESTDIR/%s' % name, objs, **kwds)
    return prog

def _AddHeaders( env, names, prefix='', stage=False ):
    if stage:
        for n in names: 
            tgt = '$INCLUDEDIR/%s/%s' % (prefix, os.path.basename(str(n)))
            env.Command( tgt, n, Copy("$TARGET", "$SOURCE") )
    return _Install(env, 'include/%s' % prefix, names) 

def _AddPythonProgram(env, name, *args, **kwds):
    wheel_prefix = 'bin'
    exes = list()
    for ver in env['PYTHONVER']:
        pyenv = env.Clone(SHOBJSUFFIX='.os%s' % ver)
        pyenv.Append(CPPFLAGS='${PYTHON%s_CPPFLAGS}' % ver)
        pylibs=env.get('PYTHON%s_LIBS' % ver)
        if pylibs: pyenv.Append(LIBS=pylibs)
        libdir = '${PYTHON%s_PREFIX}/lib' % ver
        pyenv.Append(LINKFLAGS=['-L%s' % libdir, '-Wl,-rpath,%s' % libdir])
        pyenv.Append(LINKFLAGS='${PYTHON%s_LDFLAGS}' % ver)
        objs = pyenv.AddObject(*args, **kwds)
        exe = pyenv.Program('$BINDIR/%s%s' % (name, ver), objs, **kwds) 
        _Install(env, 'bin', exe)
        exes.append(exe)
        _wheel_targets['platlib'].append((exe, wheel_prefix))
    return exes

def _AddPythonExtension(env, name, *args, **kwds):
    prefix = kwds.pop('prefix', '')
    exts = list()
    wheel_prefix = 'lib/python'
    for ver in env['PYTHONVER']:
        pyenv = env.Clone(SHOBJSUFFIX='.os%s' % ver, LDMODULESUFFIX='.so', LDMODULEPREFIX='')
        pyenv.Append(CPPFLAGS='${PYTHON%s_CPPFLAGS}' % ver)
        pylibs=env.get('PYTHON%s_LIBS' % ver)
        if pylibs: pyenv.Append(LIBS=pylibs)
        if pyenv['PLATFORM']=='darwin':
            pyenv.AppendUnique( LINKFLAGS=['-undefined','dynamic_lookup'] )
        soabi = env.get('PYTHON%s_SOABI' % ver)
        soname = name + '.%s.so' % soabi if soabi else name
        py = '$LIBDIR/python/%s/%s' % (prefix, soname)
        objs = pyenv.AddObject(*args, **kwds)
        ext = pyenv.LoadableModule(py, objs, **kwds) 
        _Install(env, 'lib/python/%s' % prefix, ext)
        exts.append(ext)
        _wheel_targets['platlib'].append((ext, wheel_prefix))
    return exts

def _AddPythonModule(env, *args, **kwds):
    exts = list()
    prefix = kwds.pop('prefix', '')
    mods = list()
    wheel_prefix = 'lib/python'
    for s in args:
        py='$LIBDIR/python/%s/%s' % (prefix, os.path.basename(s))
        mod = env.Command(py, s, Copy("$TARGET", "$SOURCE"))
        _Install(env, 'lib/python/%s' % prefix, mod)
        for ver in env['PYTHONVER']:
            cachedir = env['PYTHON%s_CACHEDIR' % ver]
            cachefile = '%s/%s/%s.%s' % (
                    os.path.dirname(py),
                    cachedir,
                    os.path.splitext(os.path.basename(py))[0],
                    env['PYTHON%s_CACHEEXT' % ver])
            cachetgt = env.Command(cachefile, mod, 'python%s -m compileall $SOURCE' % '.'.join(ver))
            _Install(env, 'lib/python/%s/%s' % (prefix, cachedir), cachetgt)
        mods.append(mod)
        _wheel_targets['purelib'].append((mod, wheel_prefix))
    return mods

def _AddShare( env, name, src=None ):
  if src is None: src=name
  if os.path.isabs(name):
      raise ValueError("AddShare: name='%s' cannot be an absolute path"%name)
  prog=env.Command('$SHAREDIR/%s'%name, src, [
      Delete("$TARGET"),
      Copy("$TARGET", "$SOURCE") ])
  _wheel_targets['data'].append((prog, 'share'))
  _InstallAs(env,'share'+'/'+name,src)
  return prog 

def _AddDoxygen( env, doxyfile ):
    prog = env.Command('$SHAREDIR/doxygen', doxyfile, '(cat $SOURCES && echo OUTPUT_DIRECTORY=$TARGET) | doxygen -')
    env.Clean(prog, prog)
    _InstallAs(env,'share/doc',prog)
    return prog

def _AddProtoC(env, src):
    objdir = env.subst('$OBJDIR')
    root = Dir('#').path
    objs = list()
    path = File(src).srcnode().get_path()
    rel = os.path.relpath(path, root)
    name = os.path.splitext(rel)[0]
    dst = '%s/%s' % (objdir, name)
    hxx, cxx =env.Command([dst+".pb.h", dst+".pb.cc"], path,
                          'protoc $SOURCE --proto_path=%s --cpp_out $OBJDIR/proto' %
                          (os.path.dirname(rel)))
    env.AddHeaders([hxx], stage=True)

    return cxx


def _AddWheel(env, tomlfile, pyver='36'):
    import enscons
    import pytoml

    with open(File(tomlfile).srcnode().abspath) as fp:
        metadata = pytoml.load(fp)['tool']['enscons']

    name = metadata['name']
    version = metadata['version']

    # obtain wheel tag using specified python version
    wmod = 'wheel' if pyver.startswith('2') else 'setuptools'
    exe = 'python%s' % '.'.join(pyver)
    tag = subprocess.check_output([exe, '-c', 'import %s.pep425tags as wp; tags=wp.get_supported(); best=[t for t in tags if "manylinux" not in "".join(t)][0]; print("-".join(best))' % wmod], universal_newlines=True).strip()

    # set things up for enscons.
    env.Replace(
            PACKAGE_NAME = name,
            PACKAGE_NAME_SAFE = name,
            PACKAGE_VERSION = version,
            PACKAGE_METADATA = metadata,
            WHEEL_TAG = tag,
            ROOT_IS_PURELIB = False,
            WHEEL_BASE = 'dist',
            DIST_BASE = 'dist',
            )
    env.Append(WHEEL_PYVER=[pyver])

    wheel_meta = enscons.init_wheel(env) 
    wheel_targets = list()
    for category, elems in _wheel_targets.items():
        if category == 'platlib':
            target_dir = env['WHEEL_PATH'].get_path()
        else:
            target_dir = env['WHEEL_DATA_PATH'].Dir(category).get_path()
        for targets, prefix in elems:
            for tgts in targets:
                for node in env.arg2nodes(tgts):
                    relpath = os.path.relpath(node.get_path(), prefix)
                    args = (os.path.join(target_dir, relpath), node)
                    wheel_targets.append(env.InstallAs(*args))

    whl = env.Zip(
            target = env['WHEEL_FILE'],
            source = wheel_meta + wheel_targets,
            ZIPROOT = env['WHEEL_PATH'])
    env.AddPostAction(whl, Action(enscons.add_manifest))
    if env.get('PREFIX'):
        out = env.Install('$PREFIX/dist/wheel', whl)
        env.Alias('install', out)
    return whl

def _Install(env,subdir,stuff,**kwds):
    if not env.get('PREFIX'): return
    tgt='$PREFIX/%s' % subdir
    if kwds.get('SHLIBVERSION'):
        EnsureSConsVersion(2,3)
        out=env.InstallVersionedLib(tgt,stuff,**kwds)
    else:
        out=env.Install(tgt,stuff,**kwds)
    env.Alias('install', tgt)
    return out

def _InstallAs(env,subdir,stuff):
    if not env.get('PREFIX'): return
    tgt = '%s/%s' % (env['PREFIX'], subdir)
    out=env.InstallAs(tgt,stuff)
    env.Alias('install', tgt)
    return out

def generate(env):
    env.Replace(
            ENV=os.environ,
            )

    opts = Variables()
    opts.Add("OBJDIR", "build product location", 'build')
    opts.Add("PREFIX", "installation location")

    opts.Add(ListVariable('PYTHONVER', 'python versions', os.getenv('PYTHONVER', ''), ['27', '35', '36', '37']))
    opts.Update(env)

    builddir = env.Dir(env['OBJDIR']).srcnode().abspath
    env['BUILDDIR'] = builddir
    env['OBJDIR'] = '$BUILDDIR/obj'
    env['LIBDIR'] = '$BUILDDIR/lib'
    env['BINDIR'] = '$BUILDDIR/bin'
    env['TESTDIR'] = '$BUILDDIR/test'
    env['SHAREDIR'] = '$BUILDDIR/share'
    env['INCLUDEDIR'] = '$BUILDDIR/include'

    # compile flags from DESRES_MODULE variables
    for key in ( 'DESRES_MODULE_CPPFLAGS',
                 'DESRES_MODULE_LDFLAGS',
                 'DESRES_MODULE_LDLIBS'):
        d=env.ParseFlags(' '.join(os.getenv(key,'').split(':')))
        env.MergeFlags(d)

    # CFLAGS and CXXFLAGS don't always get interpreted by SCons in the
    # right way (DESRESCode#1123).  Make sure the CFLAGS and CXXFLAGS
    # parts end up in the right place.
    for src, dst in (
            ('DESRES_MODULE_CFLAGS',    'CFLAGS'),
            ('DESRES_MODULE_CXXFLAGS',  'CXXFLAGS')):
        d=env.ParseFlags(' '.join(os.getenv(src,'').split(':')))
        if dst in d:
            v=d.pop(dst, None)
            if v is not None:
                env.Append(**{dst:v})
        env.MergeFlags(d)

    env.Prepend( CPPPATH=['$INCLUDEDIR'] )

    # other builders
    for name, func in globals().items():
        if name.startswith('_Add') and callable(func):
            name=name[1:]
            env.AddMethod(func, name)

    for ver in env['PYTHONVER']:
        cfg = 'python%s-config' % '.'.join(ver)
        exe = 'python%s' % '.'.join(ver)
        incs=subprocess.check_output([cfg, '--includes'], universal_newlines=True).strip()
        prefix=subprocess.check_output([cfg, '--prefix'], universal_newlines=True).strip()
        libs=subprocess.check_output([cfg, '--libs'], universal_newlines=True).strip()
        incs += ' -I' + subprocess.check_output([exe, '-c', 'import numpy; print(numpy.get_include())'], universal_newlines=True).strip()
        kwds = { 'PYTHON%s_PREFIX' % ver : prefix,
                 'PYTHON%s_CPPFLAGS' % ver : incs,
                 'PYTHON%s_LDFLAGS' % ver : libs,
                 }
        if ver.startswith('3'):
            soabi = subprocess.check_output([exe, "-c", "import sysconfig;print(sysconfig.get_config_var('SOABI'))"], universal_newlines=True).strip()
            kwds['PYTHON%s_SOABI' % ver] = soabi
            cache = subprocess.check_output([exe, '-c', 'import importlib.util as i;print(i.cache_from_source("foo.py"))'], universal_newlines=True).strip()
            kwds['PYTHON%s_CACHEDIR' % ver] = os.path.dirname(cache)
            kwds['PYTHON%s_CACHEEXT' % ver] = '.'.join(cache.split('.')[-2:])
        else:
            kwds['PYTHON%s_CACHEDIR' % ver] = '.'
            kwds['PYTHON%s_CACHEEXT' % ver] = 'pyc'

        env.Replace(**kwds)

    if env.get('PREFIX'):
        env.Alias('install', '$PREFIX')
        env.Append(RPATH=['$PREFIX/lib'], LIBPATH=['$PREFIX/lib'])
        env.Append(CPPPATH=['$PREFIX/include'])
    else:
        env.Append(RPATH=['$LIBDIR'], LIBPATH=['$LIBDIR'])

def exists(env):
    return env

def build(**kwds):
    env = Environment(**kwds)
    generate(env)
    Export('env')
    env.SConscript('SConscript', variant_dir=env['BUILDDIR'], duplicate=0)


