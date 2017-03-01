#include "system.hxx"
 
using namespace desres::msys;

int main(int argc, char *argv[]) {
    auto mol = System::create();
    auto table = mol->addTable("foo", 1);
    auto params = table->params();
    params->addParam();
    params->addProp("p", IntType);
    params->value(0,0)=42;

    return 0;
}
