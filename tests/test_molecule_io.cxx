#include "sdf.hxx"

using namespace desres::msys;

int main(int argc, char *argv[]) {
    for (int i=1; i<argc; i++) {
        auto scanner = ScanSdf(argv[i]);
        unsigned n=0;
        for (;;) {
            auto mol = scanner->next();
            if (!mol) break;
            ++n;
            auto sdf = FormatSdf(*mol);
            fwrite(sdf.data(), sdf.size(), 1, stdout);
        }
        fprintf(stderr, "%s: %u\n", argv[i], n);
    }
    return 0;
}
