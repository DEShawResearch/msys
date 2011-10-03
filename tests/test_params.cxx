#include "param_table.hxx"
#include <assert.h>
#include <stdio.h>

using namespace desres::msys;

int main(int argc, char *argv[]) {

    ParamTablePtr p = ParamTable::create();

    assert(p->paramCount()==0);
    assert(p->propCount()==0);

    /* add some properties */
    p->addProp("fc", FloatType);
    Id r0=p->addProp("r0", FloatType);
    r0=p->addProp("r0", FloatType);
    p->addProp("foo", IntType);

    assert(p->paramCount()==0);
    assert(p->propCount()==3);

    /* instantiate some params */
    Id param = p->addParam();
    assert(param==0);
    assert(p->propType(0)==FloatType);
    assert(p->propType(1)==FloatType);
    assert(p->propType(2)==IntType);
    assert(p->propName(0)=="fc");
    assert(p->propName(1)=="r0");
    assert(p->propName(2)=="foo");

    assert(p->value(param,0).asFloat()==0.0);
    assert(p->value(param,r0).asFloat()==0.0);
    assert(p->value(param,2).asFloat()==0.0);
    assert(p->value(param,2).asInt()==0.0);

    p->value(param,0).fromFloat(435.5);
    assert(p->value(param,0).asFloat()==435.5);

    p->value(param,0).fromInt(200);

    assert(p->value(param,0).asFloat()==200);
    assert(p->value(param,0)==200);
    assert(p->value(param,0)==200.);
    assert(p->value(param,0)==200.f);

    ValueRef ref = p->value(param,0);
    assert(ref==200);

    /* check that our refs aren't invalidated by additions */
    for (int i=0; i<1000; i++) p->addParam();
    assert(p->value(param,0)==200);
    assert(ref==200);

    for (int i=0; i<100; i++) {
        char buf[32];
        sprintf(buf, "%d", i);
        p->addProp(buf, StringType);
    }
    assert(p->value(param,0)==200);
    assert(ref==200);

    assert(p->value(param,0)!=201.f);

    /* test direct assignment */
    p->value(param,0)=32;
    assert(p->value(param,0)==32);

    return 0;
}
