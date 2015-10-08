#include "molfile/dtrplugin.hxx"

static void usage(FILE* fp) {
    fprintf(fp, "fsdump [--begin=n] [--end=n] [--match=xxx] [--max=n] framesetdir framesetdir...\n");
}

namespace {
    unsigned long begin=0;
    unsigned long end=-1;
    std::vector<std::string> match;
    uint64_t max=0;
    bool match_all=true;
}

static void parse_cmdline( int *pargc, char ***pargv ) {
    int i,j=0;
                           
    for (i=0; i<pargc[0]; i++) {
        pargv[0][j] = pargv[0][i];
        if ( !strcmp(pargv[0][j], "-h") ||
             !strcmp(pargv[0][j], "--help")) {
            usage(stdout);
            exit(0);
        } else if ( !strncmp(pargv[0][j], "--begin=", 8)) {
            begin = atol(pargv[0][j]+8);
        } else if ( !strncmp(pargv[0][j], "--end=", 6)) {
            end = atol(pargv[0][j]+6);
        } else if ( !strncmp(pargv[0][j], "--max=", 6)) {
            max = atol(pargv[0][j]+6);
        } else if ( !strncmp(pargv[0][j], "--match=", 8)) {
            match_all=false;
            match.push_back(pargv[0][j]+8);
        }  else {
            j++;
        }
    }
    pargc[0]=j;
    pargv[0][j]=NULL;
}


int main(int argc, char *argv[]) {
    parse_cmdline(&argc, &argv);
    void* buf = NULL;
    std::vector<int32_t> i32;
    std::vector<uint32_t> u32;
    std::vector<int64_t> i64;
    std::vector<uint64_t> u64;
    std::vector<float> f32;
    std::vector<double> f64;

    for (int i=1; i<argc; i++) {
        desres::molfile::DtrReader r(argv[i]);
        r.init();
        unsigned long size = r.size();
        unsigned long last = std::min(size, end);
        std::vector<double> times(last-begin);
        r.times(begin, times.size(), &times[0]);
        printf("{\n");
        printf("  size=%lu\n", size);
        for (unsigned long i=begin; i<last; i++) {
            printf("  time=%.17g {\n", times.at(i-begin));
            auto map = r.frame(i-begin, NULL, &buf);
            for (auto it=map.begin(), e=map.end(); it!=e; ++it) {
                if (!match_all && std::find(match.begin(), match.end(), it->first)==match.end()) continue;

                uint64_t count = it->second.count;
                printf("    %s[%lu x %s]= ", 
                        it->first.c_str(),
                        (unsigned long)count,
                        desres::molfile::dtr::Key::type_name(it->second.type));
                unsigned long nelem = max>0 ? std::min(max, it->second.count)
                                            : it->second.count;
                switch (it->second.type) {
                    case desres::molfile::dtr::Key::TYPE_INT32:
                        i32.resize(count);
                        it->second.get(&i32[0]);
                        i32.resize(nelem);
                        for (auto x : i32) printf("%d ", x);
                        break;
                    case desres::molfile::dtr::Key::TYPE_UINT32:
                        u32.resize(count);
                        it->second.get(&u32[0]);
                        u32.resize(nelem);
                        for (auto x : u32) printf("%u ", x);
                        break;
                    case desres::molfile::dtr::Key::TYPE_FLOAT32:
                        f32.resize(count);
                        it->second.get(&f32[0]);
                        f32.resize(nelem);
                        for (auto x : f32) printf("%.8g ", x);
                        break;
                    case desres::molfile::dtr::Key::TYPE_FLOAT64:
                        f64.resize(count);
                        it->second.get(&f64[0]);
                        f64.resize(nelem);
                        for (auto x : f64) printf("%.17g ", x);
                        break;
                    case desres::molfile::dtr::Key::TYPE_CHAR:
                        printf("\"%s", it->second.toString().substr(0,nelem).c_str());
                        break;

                    default:
                        ;
                }
                if (nelem < count) printf("...");
                if (it->second.type==desres::molfile::dtr::Key::TYPE_CHAR) {
                    printf("\"");
                }
                printf("\n");
            }
            printf("  }\n");
        }
        printf("}\n");
    }
    if (buf) free(buf);
    return 0;
}
