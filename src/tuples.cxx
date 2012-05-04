#include "tuples.hxx"
#include "system.hxx"
#include <boost/lexical_cast.hpp>
#include <map>

using namespace desres::msys;

#if 0
/* construct a mapping for each atom to its param as given by the
 * fan-1 table.  Throws if an atom is mapped more than once, or has
 * BadId as param.  */
static IdList find_atomtype(TermTablePtr T) {
    if (T->atomCount()!=1) {
        MSYS_FAIL("Require atomCount of 1 from " << T->name());
    }
    IdList types(T->system()->maxAtomId(), BadId);
    for (Id i=0; i<T->maxTermId(); i++) {
        if (!T->hasTerm(i)) continue;
        Id param = T->param(i);
        Id atom = T->atom(i,0);
        if (bad(param)) {
            MSYS_FAIL(T->name() << " param is undefined for atom " << atom);
        }
        if (!bad(types.at(atom))) {
            MSYS_FAIL(T->name() << " param is multiply defined for atom " << atom);
        }
        types.at(atom)=param;
    }
    return types;
}
#endif

void
desres::msys::CreateTuplesFromCombined( TermTablePtr base, 
                                        ParamTablePtr combined,
                                        TermTablePtr tuples ) {
    /* arity of tuples */
    const Id N = tuples->atomCount();
    /* columns in combined with param1, param2, ... */
    IdList pcols(N, BadId);
    /* Mapping of columns in combined to columns in base.params */
    IdList bcols(combined->propCount(), BadId);
    /* shared parameter table */
    ParamTablePtr params = base->params();

    /* map columns of combined to either param1, param2.. or base properties */
    for (Id i=0; i<combined->propCount(); i++) {
        std::string prop = combined->propName(i);
        if (prop.substr(0,5)=="param") {
            Id p = boost::lexical_cast<Id>(prop.substr(5));
            pcols.at(p-1)=i;
        } else {
            /* must be a column in base */
            Id col = params->propIndex(prop);
            if (bad(col)) {
                MSYS_FAIL("Property " << prop << " not found in " << base->name());
            }
            bcols[i] = col;
        }
    }
    if (std::find(pcols.begin(), pcols.end(), BadId)!=pcols.end()) {
        MSYS_FAIL("arity of combined table does not match " << tuples->atomCount() << " of table " << tuples->name());
    }

    /* create mapping from base type to atoms with that type */
    std::map<Id,IdList> typemap;
    for (Id i=0; i<base->maxTermId(); i++) {
        if (!base->hasTerm(i)) continue;
        Id param = base->param(i);
        typemap[param].push_back(base->atom(i,0));
    }

    /* ensure entries in combined are unique */
    std::set<IdList> combined_params;

    /* generate tuples from each entry in combined */
    for (Id i=0; i<combined->paramCount(); i++) {
        Id newp = params->addParam();
        /* copy properties from combined to base params */
        for (Id j=0; j<combined->propCount(); j++) {
            if (bad(bcols[j])) continue;
            params->value(newp,bcols[j]) = combined->value(i,j);
        }
        /* fetch param ids for this entry */
        IdList plist(N), index(N);
        Id ntuples = 1;
        for (Id j=0; j<N; j++) {
            Id param = combined->value(i,pcols[j]).asInt();
            plist[j] = param;
            ntuples *= typemap[param].size();
        }
        if (ntuples==0) {
            MSYS_FAIL("Row " << i << " in combined table references nonexistent params in " << base->name());
        }
        if (!combined_params.insert(plist).second) {
            MSYS_FAIL("Row " << i << " in combined table is a duplicate");
        }
        for (Id ituple=0; ituple<ntuples; ituple++) {
            IdList tuple(N);
            for (Id j=0; j<N; j++) tuple[j] = typemap[plist[j]][index[j]];
            tuples->addTerm(tuple,newp);
            /* advance to the next tuple */
            for (Id j=0; j<N; j++) {
                ++index[j];
                if (index[j]!=typemap[plist[j]].size()) break;
                index[j]=0;
            }
        }
    }

}

ParamTablePtr CreateCombinedFromTuples( TermTablePtr base,
                                        TermTablePtr tuples ) {

    ParamTablePtr combined = ParamTable::create();

    return combined;
}

