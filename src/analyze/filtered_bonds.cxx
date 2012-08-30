#include <boost/foreach.hpp>
#include "filtered_bonds.hxx"

namespace desres { namespace msys {

    /* remove virtuals/drudes (anything with atomic_number<1) to allow bond_order assigner and 
       other "chemistry" things to work with already typed systems */
    IdList filteredBondsForAtom(SystemPtr sys, Id aid){
        IdList filtered;
        BOOST_FOREACH(Id bid, sys->bondsForAtom(aid)){
            Id other=sys->bond(bid).other(aid);
            if(sys->atom(other).atomic_number<1) continue;
            filtered.push_back(bid);
        }
        return filtered;
    }
    IdList filteredBondedAtoms(SystemPtr sys, Id aid){
        IdList filtered;
        BOOST_FOREACH(Id bid, sys->bondsForAtom(aid)){
            Id other=sys->bond(bid).other(aid);
            if(sys->atom(other).atomic_number<1) continue;
            filtered.push_back(other);
        }
        return filtered;
    }
    Id filteredBondCountForAtom(SystemPtr sys, Id aid){
        return filteredBondsForAtom(sys,aid).size();
    }

}}
