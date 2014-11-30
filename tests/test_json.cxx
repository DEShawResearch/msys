/* @COPYRIGHT@ */

#include "fastjson/fastjson.hxx"
#include <fstream>
#include <cassert>
#include <sstream>
#include <math.h>

using desres::msys::fastjson::Json;

/* embedded quotes get written back out
 * 1.0 gets written as 1.0 and not integer 1.
 * true/false preserved 
 * keys with embedded quotes
 */

static const char js1[] = 
"[60.3, -0.001, 0., -0.0, 0.01, 0.1, -0.333, 1.0, 1.37, 60.000, true, false, { \"quote \\\"key\" : 32 }, \"embedded \\\"quotes\\\"\"]";

static void unit_test() {
    std::stringstream ss(js1);
    Json js;
    std::cout << "parsing: " << ss.str() << std::endl;
    parse_json(ss, js);

    std::stringstream ss2;
    ss2 << js << std::endl;

    std::cout << "parsing: " << ss2.str() << std::endl;
    Json js2;
    parse_json(ss2, js2);
    assert(js2.equals(js));
}

int main(int argc, char *argv[]) {
    int i;
    for (i=1; i<argc; i++) {
        std::ifstream in(argv[i]);
        Json js;
        in >> js;
        std::cout << js << std::endl;

        std::cout << "------" << std::endl;
        Json tmp;
        assert(!tmp);
        tmp.copy(js);
        std::cout << tmp << std::endl;
        assert(tmp.equals(js));
        assert(!tmp.equals(js.elem(0)));
        assert(!!tmp);
        assert(tmp.valid());
    }

    unit_test();

    return 0;
}

