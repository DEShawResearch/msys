#include "mae/mae.hxx"
#include "types.hxx"
#include <fstream>

using namespace desres::msys;
using desres::fastjson::Json;

int main(int argc, char *argv[]) {
    for (int i=1; i<argc; i++) {
        std::ifstream in(argv[i]);
        if (in) {
            double t=-now();
            Json js;
            mae::import_mae(in,js);
            t += now();
            printf("%8.3fms\n", 1000*t);
        }
    }
    return 0;
}

