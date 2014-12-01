Import('env')

# Avoid doing dependency checks on garden files.
cpp=[]
flg=[]
for p in env['CPPPATH']:
    if p.startswith('/proj'):
        flg.append('-I%s' % p)
    else:
        cpp.append(p)
env.Replace(CPPPATH=cpp)
env.Append(CFLAGS=flg, CXXFLAGS=flg)
# The install helper brings in libraries we don't want
for lib in 'molfile', 'python2.7':
    try: env.get('LIBS', []).remove(lib)
    except ValueError: pass

env.Append(
        # SSE2 for src/within.hxx.  It's optional, but way way slower without.
        CCFLAGS='-O2 -g -Wall -msse2',
        CXXFLAGS='-std=c++03 -Werror',
        )

env.SConsignFile('%s/.sconsign' % (env['OBJDIR'].strip('#')))

env.SConscript('src/SConscript')
env.SConscript('tests/SConscript')

opts=Variables()
opts.Add("MSYS_WITHOUT_PYTHON", "without python?")
opts.Update(env)

if env.get("MSYS_WITHOUT_PYTHON"):
    print "MSYS_WITHOUT_PYTHON set; will not build python extensions or tools."
else:
    env.SConscript('python/SConscript')
    env.SConscript('tools/SConscript')

