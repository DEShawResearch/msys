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

    MultiIdList fragments;
    mol->updateFragids(&fragments);

    for (unsigned frag=0; frag<fragments.size(); ++frag){
        int qTarget=FragmentChargeFromFormalCharge(mol,fragments[frag]);
        if(qTarget==INT_MAX) qTarget=0;
        if(!ChargeFromBondOrderValid(mol,fragments[frag])){
           ExportSdf( mol, "qDiff1.sdf", SdfExport::Append);
        } 
        AssignBondOrderAndFormalCharge(mol,fragments[frag]);
        ExportSdf( mol, "qAssigned.sdf", SdfExport::Append);

        IdList::iterator miter;
        printf("ATOMS:\n");
        int qTot=0;
        for (miter=fragments[frag].begin(); miter !=fragments[frag].end(); ++miter){
            atom_t const& atm1 = mol->atom(*miter);
            printf("   ( %3d %3s ): molecule= %2d  fc= %d  rc= %e  x= %f y= %f z= %f\n",
                   *miter, atm1.name.c_str(), 
                   frag, atm1.formal_charge, 
                   atm1.resonant_charge,
                   atm1.x, atm1.y, atm1.z);
            qTot+=atm1.formal_charge;
        }
        printf("BONDS:\n");
        for (miter=fragments[frag].begin(); miter !=fragments[frag].end(); ++miter){
            Id id = *miter;
            atom_t const& atm1 = mol->atom(id);
            IdList const& bonds = mol->bondsForAtom(id);
            for (unsigned j=0; j<bonds.size(); ++j){
                Id otherid = mol->bond(bonds[j]).other(id);
                if (id>otherid) continue;
                atom_t const& atm2 = mol->atom(otherid);
                double dx=atm2.x - atm1.x;
                double dy=atm2.y - atm1.y;
                double dz=atm2.z - atm1.z;
                double r=sqrt(dx*dx+dy*dy+dz*dz);
                printf("   ( %3d %s %u ) <-> ( %3d %s %u ):   order= %2.1f rorder= %f length= %f\n",
                       id,      atm1.name.c_str(), mol->bondCountForAtom(id),
                       otherid, atm2.name.c_str(), mol->bondCountForAtom(otherid),
                       (double)mol->bond(bonds[j]).order, 
                       mol->bond(bonds[j]).resonant_order,
                       r);
            }
        }
        if(qTarget!=qTot){
            printf("qDiff2: %d != %d\n",qTarget,qTot);
            ExportSdf( mol, "qDiff2.sdf", SdfExport::Append);
        }
    }
}
 
