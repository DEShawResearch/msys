
#include "system.hxx"
#include <boost/variant/get.hpp>

using namespace desres::msys;

int main() {
    auto mol = System::create();
    auto table = mol->addTable("foo", 2);
    auto& tp = table->tableProps();
    tp["nerfcut"] = 42.5;
    Float val = boost::get<Float>(tp["nerfcut"]);
    return val != 42.5;

    return 0;
}


