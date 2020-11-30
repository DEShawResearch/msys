#include "provenance.hxx"
#include <sstream>
#include <time.h>

#include <msys/version.hxx>

using namespace desres::msys;

#ifndef _MSC_VER
#include <pwd.h>
#include <unistd.h>
#include <sys/param.h>
#ifdef __APPLE__
static char* get_current_dir_name() {
    static char buf[MAXPATHLEN];
    return getwd(buf);
}
#endif
#else
static char* get_current_dir_name() {
    return ".";
}
#endif

static std::string executable_path(const char* argv0) {
#ifndef _MSC_VER
    char buf[1024];
    ssize_t rc = readlink("/proc/self/exe", buf, sizeof(buf));
    if (rc<0 || rc==1024) return argv0;
    buf[rc]='\0';
    return buf;
#else
    return argv0;
#endif
}

Provenance Provenance::fromArgs(int argc, char *argv[]) {
    Provenance prov;

    /* version */
    prov.version = "msys/";
    prov.version += MSYS_VERSION;

    /* timestamp */
    {
        char buf[200];
        time_t t;
        struct tm *tmp;
        t = time(NULL);
        tmp = localtime(&t);
        if (tmp && strftime(buf, sizeof(buf), "%c", tmp)!=0) {
            prov.timestamp = buf;
        }
    }

#ifndef _MSC_VER
    /* user */
    {
        int id = getuid();
        struct passwd* pw=getpwuid(id);
        if (pw) {
            std::stringstream ss;
            ss << id << ":" << pw->pw_name
                     << ":" << pw->pw_gecos;
            prov.user = ss.str();
        }
    }
#endif

    /* workdir */
    prov.workdir = get_current_dir_name();

    /* cmdline */
    for (int i=0; i<argc; i++) {
        prov.cmdline += argv[i];
        if (i!=argc-1) prov.cmdline += " ";
    }

    /* executable */
    if (argc > 0) {
        prov.executable = executable_path(argv[0]);
    }

    return prov;
}

const char* desres::msys::msys_version() { return MSYS_VERSION; }
int desres::msys::msys_major_version() { return MSYS_MAJOR_VERSION; }
int desres::msys::msys_minor_version() { return MSYS_MINOR_VERSION; }
int desres::msys::msys_micro_version() { return MSYS_MICRO_VERSION; }

