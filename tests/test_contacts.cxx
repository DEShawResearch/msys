#include "io.hxx"
#include "atomsel.hxx"
#include "spatial_hash.hxx"
#include "contacts.hxx"
#include <stdio.h>
#include <cmath>

using namespace desres::msys;

struct output {
    SpatialHash::ContactList& contacts;
    output(SpatialHash::ContactList& c) : contacts(c) {}
    bool exclude(Id i, Id j) const { return false; }
    void operator()(Id i, Id j, float d2) const {
        contacts.emplace_back(i,j,d2);
    }
};

int main(int argc, char *argv[]) {
    auto mol = Load(argv[1]);
    float r = atof(argv[2]);
    auto A = Atomselect(mol, argv[3]);
    auto B = Atomselect(mol, argv[4]);

    std::vector<float> pos;
    mol->getPositions(std::back_inserter(pos));
    printf("r %.3f A %lu B %lu\n", r, A.size(), B.size());
    printf("find_contacts:\n");
    double t0=now();
    SpatialHash::ContactList c1;
    find_contacts(r, &pos[0],
                  A.begin(), A.end(),
                  B.begin(), B.end(),
                  output(c1));
    double t1=now();
    printf("%u in %.3fms\n", (Id)c1.size(), (t1-t0)*1000);

    printf("spatial_hash:\n");
    double t2=now();
    SpatialHash h(&pos[0], B.size(), &B[0], NULL);
    auto c2 = h.findContacts(r, &pos[0], A.size(), &A[0]);
    double t3=now();
    printf("%u in %.3fms\n", (Id)c2.size(), (t3-t2)*1000);
    std::sort(c1.begin(), c1.end());
    std::sort(c2.begin(), c2.end());
    assert(c1==c2);

    printf("spatial_hash periodic:\n");
    double t4=now();
    SpatialHash h2(&pos[0], B.size(), &B[0], mol->global_cell[0]);
    auto c3 = h2.findContacts(r, &pos[0], A.size(), &A[0]);
    double t5=now();
    printf("%u in %.3fms\n", (Id)c3.size(), (t5-t4)*1000);
    return 0;
}

