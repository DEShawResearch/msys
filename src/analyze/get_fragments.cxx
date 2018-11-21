#include "get_fragments.hxx"
#include <map>
#include <cstdio>
#include <stack>

//#include <profiler/profiler.hxx>

namespace desres { namespace msys {

    void get_fragments(SystemPtr mol,
                       IdList const& atomIds, 
                       MultiIdList& fragments) {
        //static desres::profiler::Symbol _("desres::viparr::get_fragments"); 
        //desres::profiler::Clock __(_);

        typedef std::map<Id, Id> indexmap_t;

        fragments.clear();
        Id natoms=atomIds.size();
        if (natoms==0) return; 

        /* Create local storage for fragment Ids
         *   atoms that will have fragment Ids assigned == BadId */
        IdList assignments(atomIds.size(),BadId);
        indexmap_t atom_to_index;
        indexmap_t::iterator iter=atom_to_index.end();
        for (Id i=0, n=atomIds.size(); i<n; ++i){
            Id aid=atomIds[i];
            assert(mol->hasAtom(aid));
            iter=atom_to_index.insert(iter,
                    indexmap_t::value_type(aid, i));
        }

        std::stack<Id> S;
        Id fragid=0;
        for (Id i=0, n=atomIds.size(); i<n; i++) {
            if (assignments[i]!=BadId) continue;   /* already assigned */
            S.push(atomIds[i]);
            assignments[i] = fragid;
            do {
                Id aid=S.top();
                S.pop();
                const IdList& bonds = mol->bondsForAtom(aid);
                for (IdList::const_iterator j=bonds.begin(), e=bonds.end(); j!=e; ++j) {
                    Id other = mol->bond(*j).other(aid);
                    iter = atom_to_index.find(other);
                    if (iter != atom_to_index.end() && assignments[iter->second]==BadId){
                        assignments[iter->second]=fragid;
                        /* Only add this atom if its non-terminal */
                        if(mol->bondCountForAtom(other)>1) S.push(other);
                    }
                }
            } while (S.size());
            ++fragid;
        }
        
        fragments.resize(fragid);
        for (Id idx = 0; idx < natoms; ++idx){
            fragments[assignments[idx]].push_back(atomIds[idx]);
        }
    }

}}

    
