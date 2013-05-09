#include <stdio.h>
#include <math.h>
#include "analyze.hxx"
#include "clone.hxx"
//#include <profiler/profiler.hxx>
#include "load.hxx"
#include "sdf.hxx"
#include "atomsel.hxx"

using namespace desres::msys;

static void assign(SystemPtr);

int main(int argc,char **argv){

    for (int i=1; i<argc; i++) {
            LoadIteratorPtr it=LoadIterator::create(argv[i],true);
            SystemPtr mol;
            Id count=0;
            while (mol=it->next()){
               mol = Clone(mol, Atomselect(mol, "atomicnumber > 0"));
            try {
               printf("Entry %u\n",count);
               assign(mol);
            } catch (std::exception& e) {
               fprintf(stderr, "  FAILED ON %s: %s\n", argv[i], e.what());
               ExportSdf( mol, "qFail.sdf", SdfExport::Append);
            }
            count++;
        }
    }
    //desres::profiler::printFlat(stdout);
    //desres::profiler::printTree(stdout);
    return 0;
}

static
bool ChargeFromBondOrderValid(SystemPtr mol,
                              IdList const& atoms){

    int qSum=FragmentChargeFromAtomicCharge(mol, atoms);
    int fcSum=FragmentChargeFromFormalCharge(mol, atoms);
    int bqSum=FragmentChargeFromBondOrders(mol, atoms);

    bool good=false;
    if(fcSum!=INT_MAX){
        good=(fcSum==bqSum);
    }else if(qSum!=INT_MAX){
        good=(qSum==bqSum);
    }else if(bqSum!=INT_MAX){
        good=true;
    }

    if(good){
        printf("qSummary+: %d %d %d\n",qSum, fcSum, bqSum);
    }else{
        printf("qSummary-: %d %d %d\n",qSum, fcSum, bqSum);
    }
    printf("Total charge detected via ff atom charge: %d\n",qSum);
    printf("Total charge detected via formal charge : %d\n",fcSum);
    printf("Total charge detected via bondOrders    : %d\n",bqSum);

    return good;
}


static
void assign(SystemPtr mol) {

    IdList topids=ComputeTopologicalIds(mol);
    MultiIdList fragments;
    mol->updateFragids(&fragments);

    Id badQ=0;
    int qTarget=0;
    int qTot=0;
    std::vector<int> qTargets;
    std::map<Id, std::map<int, int> > tid_to_qCount;
    BOOST_FOREACH(Id aid, mol->atoms()){
        int q=mol->atom(aid).formal_charge;
        qTarget+=q;
        qTargets.push_back(q);
        tid_to_qCount[topids[aid]][q]+=1;
    }

    try{
        for (unsigned frag=0; frag<fragments.size(); ++frag){
            int qTarget=FragmentChargeFromFormalCharge(mol,fragments[frag]);
            if(qTarget==INT_MAX) qTarget=0;
            if(false && !ChargeFromBondOrderValid(mol,fragments[frag])){
                ExportSdf( mol, "qDiff1.sdf", SdfExport::Append);
            } 
            
            AssignBondOrderAndFormalCharge(mol,fragments[frag]);
            ExportSdf( mol, "qAssigned.sdf", SdfExport::Append);
            
            
            BOOST_FOREACH(Id aid, fragments[frag]){
                int q=mol->atom(aid).formal_charge;
                std::map<int, int> &qTargets=tid_to_qCount.find(topids[aid])->second;
                std::map<int, int>::iterator iter=qTargets.find(q);
                if(iter==qTargets.end() || iter->second==0){
                    badQ++;
                }else{
                    iter->second--;
                }
                qTot+=q;
            }   
        }
        if(badQ){
            ExportSdf( mol, "qAtomDiff.sdf", SdfExport::Append);
        }
        if(qTot != qTarget){
            ExportSdf( mol, "qTotDiff.sdf", SdfExport::Append);
        }
    } catch (std::exception& e) {
        fprintf(stderr, "  FAILED %s\n", e.what());
        ExportSdf( mol, "qFail.sdf", SdfExport::Append);
    }

#if 1
    for (unsigned frag=0; frag<fragments.size(); ++frag){
        printf("Frag %u ATOMS:\n",frag);
        BOOST_FOREACH(Id aid, fragments[frag]){
            atom_t const& atm1 = mol->atom(aid);
            printf("   ( %3d %3s ): molecule= %2u  fc= %d  rc= %e  x= %f y= %f z= %f\n",
                   aid, atm1.name.c_str(), 
                   frag, atm1.formal_charge, 
                   atm1.resonant_charge,
                   atm1.x, atm1.y, atm1.z);
        }

        printf("Frag %u BONDS:\n",frag);
        BOOST_FOREACH(Id aid, fragments[frag]){
            atom_t const& atm1 = mol->atom(aid);
            IdList const& bonds = mol->bondsForAtom(aid);
            for (unsigned j=0; j<bonds.size(); ++j){
                Id otherid = mol->bond(bonds[j]).other(aid);
                if (aid>otherid) continue;
                atom_t const& atm2 = mol->atom(otherid);
                double dx=atm2.x - atm1.x;
                double dy=atm2.y - atm1.y;
                double dz=atm2.z - atm1.z;
                double r=sqrt(dx*dx+dy*dy+dz*dz);
                printf("   ( %3d %s %u ) <-> ( %3d %s %u ):   order= %2.1f rorder= %f length= %f\n",
                       aid,      atm1.name.c_str(), mol->bondCountForAtom(aid),
                       otherid, atm2.name.c_str(), mol->bondCountForAtom(otherid),
                       (double)mol->bond(bonds[j]).order, 
                       mol->bond(bonds[j]).resonant_order,
                       r);
            }
        }
    }
#endif
}
 
