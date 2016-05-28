#include "io.hxx"
#include "atomsel.hxx"
#include "spatial_hash.hxx"
#include "contacts.hxx"
#include <stdio.h>

using namespace desres::msys;

struct output {
    bool exclude(Id i, Id j) const { return false; }
    void operator()(Id i, Id j, float d2) const {
        printf("%u %u %f\n", i,j,sqrt(d2)); 
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
    find_contacts(r, &pos[0],
                  A.begin(), A.end(),
                  B.begin(), B.end(),
                  output());
    double t1=now();
    printf("%.3fms\n", (t1-t0)*1000);

    printf("spatial_hash:\n");
    double t2=now();
    SpatialHash h(&pos[0], B.size(), &B[0], NULL);
    for (auto c : h.findContacts(r, &pos[0], A.size(), &A[0])) {
        printf("%u %u %f\n", c.i, c.j, c.d);
    }
    double t3=now();
    printf("%.3fms\n", (t3-t2)*1000);

    printf("spatial_hash periodic:\n");
    double t4=now();
    SpatialHash h2(&pos[0], B.size(), &B[0], mol->global_cell[0]);
    for (auto c : h2.findContacts(r, &pos[0], A.size(), &A[0])) {
        printf("%u %u %f\n", c.i, c.j, c.d);
    }
    double t5=now();
    printf("%.3fms\n", (t5-t4)*1000);
    return 0;
}

