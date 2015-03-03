#include <stdint.h>
#include <deque>
#include <sstream>
#include <cstdio>

#include <boost/foreach.hpp>

#include "../param_table.hxx"
#include "../analyze.hxx"
#include "../graph.hxx"

using namespace desres::msys;

typedef std::pair<uint32_t,uint32_t> tid_t;
typedef std::map<tid_t, std::vector<Id> > tidmap_t;

static const Id DEFAULTTID = 0;

namespace {
    uint32_t atomic_invariant(const SystemPtr mol, Id aid1){

        IdList const& bonds = mol->bondsForAtom(aid1);
        uint32_t nbonds=bonds.size();
        uint32_t nterminal=0;
        uint32_t aenv=1;
        BOOST_FOREACH(Id bid, bonds){
            Id aid2= mol->bond(bid).other(aid1);
            if(mol->bondCountForAtom(aid2)==1) nterminal+=1;
            int anum=mol->atom(aid2).atomic_number;
            switch (anum){
            case 1:  // H
                break;
            case 6:  // C
                aenv*=2;
                break;
            case 7:  // N
                aenv*=3;
                break;
            case 8:  // O
                aenv*=5;
                break;
            case 15: // P
                aenv*=7;
                break;
            case 16: // S
                aenv*=11;
                break;
            default: // everything else
                aenv*=13;             
                break;
            }
        }
        uint32_t atomprop =            // MAXBONDS = 5 anticipated
            (nbonds           << 0)  | // 3  bits (MAXBONDS < 8)
            (nterminal        << 3)  | // 3  bits (MAXBONDS < 8)
            (mol->atom(aid1).atomic_number << 6)  | // 7  bits (atomic # <= 127 )
            (aenv             << 13) ; // 19 bits ( 13**MAXBONDS  = 371293 < 524288 [2**19] )
        return atomprop;
    }
    
    uint32_t update_tids_from_tidmap(tidmap_t const& fragTids,
                                     std::vector<uint32_t> const& primes,
                                     std::vector<uint32_t> &alltids){
        typedef std::map<tid_t, uint32_t> cache_t;
        
        if(fragTids.size()==0 ) return 0;
        cache_t cache;
        
        uint32_t ntid=DEFAULTTID;
        for(tidmap_t::const_iterator titer=fragTids.begin(); titer != fragTids.end(); ++titer){
            cache_t::iterator citer=cache.lower_bound(titer->first);
            if(citer == cache.end() || cache.key_comp()(titer->first,citer->first)){
                ntid+=1;
                citer=cache.insert(citer, cache_t::value_type(titer->first,primes.at(ntid)));
            }
            BOOST_FOREACH(Id aid1, titer->second){
                alltids[aid1]=citer->second;
            }
        }
        return ntid;
    }
    uint32_t finalize_tids(SystemPtr mol, IdList const& atomIds,
                           tidmap_t const& fragTids, 
                           std::vector<uint32_t> &alltids,
                           Id tid_offset){
        
        typedef std::map<tid_t, Id> cache_t;
        
        /* generate consecutive ids to classified atoms */
        Id ntid=tid_offset;
        for(tidmap_t::const_iterator titer=fragTids.begin(); titer != fragTids.end(); ++titer){
            ntid+=1;
            BOOST_FOREACH(Id aid1, titer->second){
                alltids[aid1]=ntid;
            }
        }
        
        /* assign tids to unclassified (terminal) atoms */
        cache_t cache;
        BOOST_FOREACH(Id aid1, atomIds){
            if(alltids[aid1]==DEFAULTTID){
                assert(mol->bondCountForAtom(aid1)==1);
                Id aid2=mol->bondedAtoms(aid1)[0];
                tid_t tid(alltids[aid2],static_cast<uint32_t>(mol->atom(aid1).atomic_number));
                cache_t::iterator citer=cache.lower_bound(tid);
                if(citer == cache.end() || cache.key_comp()(tid,citer->first)){
                    ntid+=1;
                    citer=cache.insert(citer, cache_t::value_type(tid,ntid));
                }
                alltids[aid1]=citer->second;
            }
        }
        return ntid-tid_offset;
    }

    /* skip terminal atoms */
    inline bool skip_this_atom(const SystemPtr mol, Id aid) {
        return mol->bondCountForAtom(aid)==1; 
    }
    
    void generate_primes(Id n, std::vector<Id> &primes) {
        /* This is about the simplest prime generation algorithm I could
         * think of.  It will be called with n equal to the number of atoms
         * in the largest fragment in a system.  The largest such fragment
         * I've been able to find has 15k atoms, for which this simple sieve 
         * procedure takes around 5ms.  For n=1e6, which I would consider
         * to be highly unlikely, it takes around 1.4s, which is acceptable.
         * JRG - 3 March 2015 */
        primes.clear();
        primes.reserve(n+1);
        primes.push_back(DEFAULTTID);
        if (n>0) primes.push_back(2);
        Id candidate = 3;
        while (primes.size() <= n) {
            bool prime = true;
            for (Id k=2, m=primes.size(); k<m; k++) {
                Id p = primes[k];
                if (p*p > candidate) break;
                if (candidate % p == 0) {
                    prime = false;
                    break;
                }
            }
            if (prime) {
                primes.push_back(candidate);
            }
            candidate += 2;
        }
    }    
}

Id assign_fragment_topological_ids(SystemPtr mol, 
                                   IdList const& fragment,
                                   std::vector<uint32_t> const& primes,
                                   std::vector<uint32_t> &alltids,
                                   Id tid_offset){

    assert(primes.size()>=fragment.size()+1);
    /* tidmaps are used to map unique (class,property) tuples into topological ids 
       fragTids contains tids for the current fragment  */
    tidmap_t fragTids;
    
    /* Initalize atom tids -
       skip terminal atoms, rest are active (could freeze ions but its not necessary) */
    BOOST_FOREACH(Id aid1, fragment){
        if(skip_this_atom(mol,aid1)) continue;
        uint32_t prop = atomic_invariant(mol,aid1);
        /* Initial atom classes are (0, atom invarient) */
        tid_t tid(0,prop);
        fragTids[tid].push_back(aid1);           
    }
 
    /* The iteration count of the algorithm is bounded by the fragment size */
    Id maxloop=fragment.size();

    /* run Simple partitioning + extended partitioning until # of tids stop changing */
    Id ntid_extended=0;
    Id ntid_simple=update_tids_from_tidmap(fragTids, primes, alltids);
    assert(ntid_simple==fragTids.size());
    do {
        /* map tuples into unique tids (prime #'s) */
        for ( Id iloop=0; iloop<maxloop; ++iloop){
            /* reset our notion of fragTids */
            fragTids.clear();
        
            /* calculate new atom properties based local environment */
            BOOST_FOREACH(Id aid1, fragment){
                /* skip non branched atoms... We can add them back later */
                if(skip_this_atom(mol,aid1)) continue;
            
                Id prop=1;  
                IdList const& bonds=mol->bondsForAtom(aid1);
                BOOST_FOREACH(Id bid, bonds){
                    Id aid2 = mol->bond(bid).other(aid1);
                    if(skip_this_atom(mol, aid2)) continue;
                    prop*=alltids[aid2];
                }
                
                /* New atom classes (current tid, new atom prop)
                   This prevents the tid ping-pong effect inherent in the
                   orignal version of the morgan algorithm. */
                tid_t tid(alltids[aid1], prop);
                fragTids[tid].push_back(aid1);  
            }
            
            Id ntid_new=update_tids_from_tidmap(fragTids, primes, alltids);
            assert(ntid_new==fragTids.size());

            /* If the number of generated classes is the same as the previous iteration
               the tid assignment for this fragment is complete, otherwise its not. */
            if(ntid_simple == ntid_new){
                break;
            }
            ntid_simple = ntid_new;
        }
        
        /* short circuit for second iteration */
        if(ntid_simple==ntid_extended)break;

        /* Simple stage completed, now brute force check equivalent tids to 
           seperate out problematic cases */
        GraphPtr fGraph=Graph::create(mol, fragment);
        IdList const& gids=fGraph->atoms();
        std::vector<int> attr;
        BOOST_FOREACH(Id aid1, gids){
            if(alltids[aid1]==DEFAULTTID){
                attr.push_back(primes[ntid_simple]+mol->atom(aid1).atomic_number);
            }else{
                attr.push_back(alltids[aid1]);
            }
        }
        fGraph->setNodeAttributes(attr);
        tidmap_t extTids;
        for(tidmap_t::const_iterator titer=fragTids.begin(); titer != fragTids.end(); ++titer){
            IdList primaryNodes;
            BOOST_FOREACH(Id aid1, titer->second){
                bool matched=false;
                std::vector<IdPair> matches;
                for(uint32_t pidx=0; pidx<primaryNodes.size();++pidx){
                    if(fGraph->match(fGraph, primaryNodes[pidx], aid1, matches)){
                        matched=true;
                        tid_t tid(alltids[aid1], pidx);
                        extTids[tid].push_back(aid1);  
                        break;
                    }
                }
                if(!matched){
                    tid_t tid(alltids[aid1], primaryNodes.size());
                    extTids[tid].push_back(aid1);   
                    primaryNodes.push_back(aid1);
                }
            }
        }
        std::swap(fragTids,extTids);
        ntid_extended=update_tids_from_tidmap(fragTids, primes, alltids);
        assert(ntid_extended==fragTids.size());

    }while(ntid_simple!=ntid_extended);

    /* Finalize the atom tids */
    return finalize_tids(mol,fragment,fragTids,alltids,tid_offset);
    
}


namespace desres { namespace msys {
    /* 
       based partially on:
       SMILES. 2. Algorithm for Generation of Unique SMILES Notation
       J. Chem. Inf Comput. Sci., Vol. 29, No. 2, 1989
    */
    IdList ComputeTopologicalIds(SystemPtr mol) {
        Id maxaid=mol->maxAtomId();
        if(maxaid==0) return IdList();

        /* 0) get fragments */
        MultiIdList fragments;
        Id maxfrags=mol->updateFragids(&fragments);

        /* 1) Map fragments into possibly equivalent groups based on graph hash */
        typedef std::map<std::string, IdList> fragmap_t;
        fragmap_t fragmap;
        Id maxFrag=0;
        for(Id ifrag=0; ifrag<maxfrags;++ifrag){
            std::string hash=Graph::hash(mol,fragments[ifrag]);
            fragmap[hash].push_back(ifrag);
            if(fragments[ifrag].size() > maxFrag) maxFrag=fragments[ifrag].size();
        }

        /* 2) Initialize primes and tid storage */
        IdList primes;
        generate_primes(maxFrag, primes);
        IdList alltids(maxaid,DEFAULTTID);

        /* 3) Assign tids to uniqe fragment within a hashed group, clone tids to rest */
        Id offset=DEFAULTTID;
        for(fragmap_t::iterator frags=fragmap.begin(); frags!= fragmap.end(); ++frags){
            IdList const& fragIdlist=frags->second;
   
            std::vector<GraphPtr> primaryGraphs;
            BOOST_FOREACH(Id fragid, fragIdlist){
                bool matched=false;
                GraphPtr gCurrent=Graph::create(mol, fragments[fragid]);
                std::vector<IdPair> matches;
                BOOST_FOREACH(GraphPtr const& gPrimary, primaryGraphs){
                    if( gCurrent->match(gPrimary, matches)){
                        matched=true;
                        /* Clone tids from matched to current */
                        BOOST_FOREACH(IdPair const& p, matches){
                            alltids[p.first]=alltids[p.second];
                        }
                        break;
                    }
                }
                if(!matched){
                    /* Save graph and assign topological ids */
                    primaryGraphs.push_back(gCurrent);
                    offset+=assign_fragment_topological_ids(mol, fragments[fragid], 
                                                            primes, alltids, offset); 
                }
            }
        }

        return alltids;
    }

}}
