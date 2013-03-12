#include <dms/dms.hxx>
#include <stdio.h>

using namespace desres::msys;

int main(int argc, char *argv[]) {
    for (int i=1; i<argc; i++) {
        const char* path = argv[i];
        Sqlite dms = Sqlite::read(path);
        Reader r = dms.fetch("dms_version");
        if (!r.size()) {
            printf("%s unknown\n", path);
            continue;
        }
        int MAJOR = r.column("major");
        int MINOR = r.column("minor");
        if (MAJOR<0 || MINOR<0) {
            MSYS_FAIL("dms_version table is malformatted in " << path);
        }
        int major = r.get_int(MAJOR);
        int minor = r.get_int(MINOR);
        printf("%s %d.%d\n", path, major, minor);
    }
    return 0;
}
