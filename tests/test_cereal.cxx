#include <msys/io.hxx>
#include <msys/cereal.hxx>
#include <msys/hash.hxx>
#include <cassert>
#include <fstream>

using namespace desres::msys;

int main(int argc, char *argv[]) {
    for (int i=1; i<argc; i++) {
        FileFormat fmt;
        double t=-now();
        SystemPtr m = Load(argv[i], false, &fmt);
        t+=now();
        printf("%s %.3fs\n", argv[i], t);
        if (!m) MSYS_FAIL("Unable to guess filetype of " << argv[i]);

        std::stringstream ss;
        t=-now();
        ExportCereal(m, ss, Provenance::fromArgs(argc, argv));
        t+=now();
        printf("  export cereal %.3fs\n", t);
        t=-now();
        auto m2 = ImportCereal(ss);
        t+=now();
        printf("  import cereal %.3fs\n", t);

        std::ofstream out("out.cereal");
        ExportCereal(m, out, Provenance::fromArgs(argc, argv));
        out.close();

        assert(HashSystem(m) == HashSystem(m2));
    }
    return 0;
}




