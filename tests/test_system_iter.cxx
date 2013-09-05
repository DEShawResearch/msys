#include "dms.hxx"
#include <cstdio>
#include <boost/foreach.hpp>

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
        BOOST_FOREACH(Id id, std::make_pair(mol->atomBegin(), mol->atomEnd()))
            ids.push_back(id);
        assert(ids.at(0)==1);
        assert(ids.at(1)==3);
        assert(ids.size()==2);
    }

    return 0;
}

