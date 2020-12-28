#include "dms.hxx"
#include <cstdio>
#include <cassert>

using namespace desres::msys;

int main(int argc, char *argv[]) {
    SystemPtr mol = System::create();
    Id res=mol->addResidue(mol->addChain());
    mol->addAtom(res);
    mol->addAtom(res);
    mol->addAtom(res);
    mol->addAtom(res);

    {
        IdList ids;
        std::copy(mol->atomBegin(), mol->atomEnd(), std::back_inserter(ids));
        for (Id i=0; i<4; i++) assert(ids.at(i)==i);
    }

    mol->delAtom(2);
    {
        IdList ids;
        std::copy(mol->atomBegin(), mol->atomEnd(), std::back_inserter(ids));
        assert(ids.at(0)==0);
        assert(ids.at(1)==1);
        assert(ids.at(2)==3);
        assert(ids.size()==3);
    }

    mol->delAtom(0);
    {
        IdList ids;
        std::copy(mol->atomBegin(), mol->atomEnd(), std::back_inserter(ids));
        assert(ids.at(0)==1);
        assert(ids.at(1)==3);
        assert(ids.size()==2);
    }

    {
        IdList ids;
        for (auto id : mol->atoms()) ids.push_back(id);
        assert(ids.at(0)==1);
        assert(ids.at(1)==3);
        assert(ids.size()==2);
    }

    mol = System::create();
    res=mol->addResidue(mol->addChain());
    for (Id i=0; i<100000; i++) mol->addAtom(res);
    double t=-now();
    std::count(mol->atomBegin(), mol->atomEnd(), 0);
    t+=now();
    printf("iter: %8.3f\n", t*1000);
    t=-now();
    int n=0;
    for (Id i=0; i<mol->maxAtomId(); i++) { 
        if (!mol->hasAtom(i)) continue;
        if (i==0) ++n;
    }
    t+=now();
    printf("old: %8.3f\n", t*1000);

    return 0;
}

