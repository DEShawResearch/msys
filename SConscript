Import('env')
env.Append(
        CFLAGS='-O2 -g',
        CXXFLAGS='-O2 -Wall -Werror -g -std=c++0x',
        )

env.SConsignFile('%s/.sconsign' % (env['OBJDIR'].strip('#')))

env.SConscript('src/SConscript')
env.SConscript('tests/SConscript')

opts=Variables()
opts.Add("WITHOUT_PYTHON", "without python?")
opts.Update(env)

if env.get("WITHOUT_PYTHON"):
    print "WITHOUT_PYTHON set; will not build python extensions or tools."
else:
    env.SConscript('python/SConscript')
    env.SConscript('tools/SConscript')

