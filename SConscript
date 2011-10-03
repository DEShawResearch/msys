Import('env')
env.Append(
        CFLAGS='-O2 -g',
        CXXFLAGS='-O2 -Wall -Werror -g -std=c++0x',
        )

env.SConsignFile('%s/.sconsign' % (env['OBJDIR'].strip('#')))

env.SConscript('src/SConscript')
env.SConscript('tests/SConscript')
env.SConscript('python/SConscript')
env.SConscript('tools/SConscript')
