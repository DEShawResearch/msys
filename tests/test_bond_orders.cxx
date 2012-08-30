#include <stdio.h>
#include <math.h>
#include "analyze.hxx"
#include "clone.hxx"
//#include <profiler/profiler.hxx>
#include "dms.hxx"
#include "atomsel.hxx"

using namespace desres::msys;

int main(int argc,char **argv){

    if (argc!=3) {
        fprintf(stderr, "Usage: %s input.dms selection\n", argv[0]);
        return 2;
    }
    const char* ipath=argv[1];
    std::string sel("atomicnumber>0 and (");
    sel += argv[2];
    sel += ")";

    const bool structure_only = true;
    SystemPtr mol = ImportDMS(ipath, structure_only);
    mol = Clone(mol, Atomselect(mol, sel));

    AssignBondOrderAndFormalCharge(mol);

    MultiIdList fragments;
    mol->updateFragids(&fragments);

    for (unsigned frag=0; frag<fragments.size(); ++frag){
        IdList::iterator miter;
        printf("ATOMS:\n");
        for (miter=fragments[frag].begin(); miter !=fragments[frag].end(); ++miter){
            atom_t const& atm1 = mol->atom(*miter);
            printf("   ( %3d %3s ): molecule= %2d  fc= %d  rc= %e  x= %f y= %f z= %f\n",
                   *miter, atm1.name.c_str(), 
                   frag, atm1.formal_charge, 
                   mol->atomPropValue(*miter, "resonant_charge").asFloat(),
                   atm1.x, atm1.y, atm1.z);
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
                       mol->bond(bonds[j]).order, 
                       mol->bondPropValue(bonds[j],"resonant_order").asFloat(),
                       r);
            }
        }
    }

    //desres::profiler::printFlat(stdout);
    //desres::profiler::printTree(stdout);

    return 0; // exit with success
}
 
