#include "override.hxx"
#include "system.hxx"

using namespace desres::msys;

OverrideMap
desres::msys::FindOverridesFromTuples( TermTablePtr base,
                                       TermTablePtr tuples ) {

    /* input checks */
    if (base->atomCount()!=1) MSYS_FAIL("base must have atomCount==1");
    if (tuples->atomCount()!=2) MSYS_FAIL("tuples must have atomCount==2");

    typedef std::map<std::pair<Id,Id>, Id> OverrideMap;
    OverrideMap map;

    /* lookup from atom to parameter */
    IdList nbtype(base->system()->maxAtomId(), BadId);
    for (Id i=0; i<base->maxTermId(); i++) {
        if (!base->hasTerm(i)) continue;
        nbtype.at(base->atom(i,0)) = base->param(i);
    }

    /* iterate over tuples */
    std::pair<Id,Id> types;
    for (Id i=0; i<tuples->maxTermId(); i++) {
        if (!tuples->hasTerm(i)) continue;
        Id ai = tuples->atom(i,0);
        Id aj = tuples->atom(i,1);
        Id param = tuples->param(i);
        types.first  = nbtype.at(ai);
        types.second = nbtype.at(aj);
        if (bad(types.first))  MSYS_FAIL("no nonbonded type for atom " << ai);
        if (bad(types.second)) MSYS_FAIL("no nonbonded type for atom " << aj);
        if (types.first > types.second) std::swap(types.first, types.second);
        if (map.count(types)) {
            if (tuples->params()->compare(param, map[types])) {
                MSYS_FAIL("conflicting override for particles " << ai << ", " << aj);
            }
        } else {
            map[types] = param;
        }
    }
    return map;
}

void 
desres::msys::MakeTuplesFromOverrides( OverrideMap const& o,
                                       TermTablePtr base,
                                       TermTablePtr tuples ) {

    for (OverrideMap::const_iterator it=o.begin(); it!=o.end(); ++it) {
        const Id param1 = it->first.first;
        const Id param2 = it->first.second;
        const Id p = it->second;
        IdList ids(2);
        for (Id ti=0; ti<base->maxTermId(); ti++) {
            if (!base->hasTerm(ti)) continue;
            if (base->param(ti)!=param1) continue;
            ids[0] = base->atom(ti,0);
            if (param1==param2) {
                ids[1] = ids[0];
                tuples->addTerm(ids,p);
            } else for (Id tj=0; tj<base->maxTermId(); tj++) {
                if (!base->hasTerm(tj)) continue;
                if (base->param(tj)!=param2) continue;
                ids[1] = base->atom(tj,0);
                tuples->addTerm(ids,p);
            }
        }
    }
}

