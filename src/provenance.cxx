#include "provenance.hxx"
#include <boost/filesystem.hpp>
#include <sstream>
#include <time.h>
#include <pwd.h>

using namespace desres::msys;
namespace bfs = boost::filesystem;

Provenance Provenance::fromArgs(int argc, char *argv[]) {
    Provenance prov;

    /* timestamp */
    {
        char buf[200];
        time_t t;
        struct tm *tmp;
        t = time(NULL);
        printf("time: %f\n", (double)t);
        tmp = localtime(&t);
        printf("tmp: %p\n", (void *)tmp);
        if (tmp && strftime(buf, sizeof(buf), "%c", tmp)!=0) {
            prov.timestamp = buf;
        }
    }

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

    /* workdir */
    prov.workdir = bfs::system_complete(".").string();

    /* cmdline */
    for (int i=0; i<argc; i++) {
        prov.cmdline += argv[i];
        if (i!=argc-1) prov.cmdline += " ";
    }

    return prov;
}
