#include "value.hxx"
#include <string.h>
#include <assert.h>
#include <cstdlib>

using namespace desres::msys;

int main() {

    Int i1=1;
    Float f1=1;
    char* s1 = strdup("1");

    Int i2=2;
    Float f2=1.5;
    char* s2 = strdup("1.5");

    Value vi1, vf1, vs1;
    vi1.i=i1;
    vf1.f=f1;
    vs1.s=s1;

    Value vi2, vf2, vs2;
    vi2.i=i2;
    vf2.f=f2;
    vs2.s=s2;

    ValueRef ri1(IntType, vi1);
    ValueRef rf1(FloatType, vf1);
    ValueRef rs1(StringType, vs1);

    ValueRef ri2(IntType, vi2);
    ValueRef rf2(FloatType, vf2);
    ValueRef rs2(StringType, vs2);


    /* self equivalence */
    assert(ri1==ri1);
    assert(rf1==rf1);
    assert(rs1==rs1);

    /* comparison */
    assert(ri1!=ri2);
    assert(rf1!=rf2);
    assert(rs1!=rs2);

    /* ints and floats interconvert */
    assert(rf1==ri1);
    assert(ri1==rf1);

    assert(ri1!=rf2);
    assert(rf1!=ri2);
    assert(rf2!=ri1);
    assert(ri2!=rf1);

    /* strings don't interconvert */
    assert(rs1!=ri1);
    assert(ri1!=rs1);
    assert(rs1!=rf1);
    assert(rf1!=rs1);

    free(s1);
    free(s2);
    return 0;
}
