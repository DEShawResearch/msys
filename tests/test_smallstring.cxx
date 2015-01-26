#include "smallstring.hxx"
#include "value.hxx"

using namespace desres::msys;

int main(int argc, char *argv[]) {

    SmallString<10> x;
    std::cout <<  "something " << x << " more things";
    try {
        MSYS_FAIL("something " << x << " more things");
    } catch (...) {
    }
    return 0;
}



