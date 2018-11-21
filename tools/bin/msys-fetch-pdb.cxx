#include "pdb.hxx"
#include <stdio.h>

static void print_usage(FILE* fp, const char* argv0) {
    fprintf(fp, "Usage: msys-fetch-pdb CODE [-o ofile] [-h,--help]\n");
}

static const char* ofile = nullptr;

static std::vector<char*> parse_cmdline(int argc, char **argv) {
    std::vector<char *> args;
    for (int i=1; i<argc; i++) {
        if ( !strcmp(argv[i], "-h") ||
             !strcmp(argv[i], "--help")) {
            print_usage(stdout, argv[0]);
            exit(0);
        } else if (!strcmp(argv[i], "-o") ||
                   !strcmp(argv[i], "--output")) {
            if (i+1==argc) {
                print_usage(stderr, argv[0]);
                exit(1);
            }
            ofile = argv[++i];
        }  else {
            args.push_back(argv[i]);
        }
    }
    return args;
}

int main(int argc, char *argv[]) {
    auto args = parse_cmdline(argc, argv);
    if (args.empty()) {
        print_usage(stderr, argv[0]);
        exit(1);
    }
    FILE* fp= ofile ? fopen(ofile, "w") : stdout;
    if (!fp) {
        fprintf(stderr, "%s: Unable to open output file %s for writing\n",
                argv[0], ofile ? ofile : "<stdout>");
        exit(1);
    }
    std::shared_ptr<FILE> dtor(fp, fclose);
    for (auto arg : args) {
        auto s = desres::msys::FetchPDB(arg);
        if (fwrite(s.data(), s.size(), 1, fp) != 1) {
            fprintf(stderr, "%s: %s\n", argv[0], strerror(errno));
            exit(1);
        }
    }
    dtor.reset();
    return 0;
}


