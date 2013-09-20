#include "contacts.hxx"
#include "io.hxx"

using namespace desres::msys;

struct Output {
    bool exclude(Id i, Id j) const { return false; }
    void operator()(Id i, Id j, double r2) const {
        printf("%u %u %8.3f\n", i, j, sqrt(r2));
    }
};

int main(int argc, char *argv[]) {
    for (int arg=1; arg<argc; arg++) {
        const char* path = argv[arg];
        printf("%s\n", path);
        SystemPtr mol = Load(path);
        Id i,n = mol->atomCount();
        if (!n) continue;
        std::vector<double> pos(3*n);
        IdList A, B;
        for (i=0; i<3; i++) {
            //A.push_back(i);
            //B.push_back(i+3);
        }
        A.push_back(3030);
        B.push_back(4679);
        for (i=0; i<n; i++) {
            memcpy(&pos[3*i], &mol->atom(i).x, 3*sizeof(double));
        }
        printf("  NON-PERIODIC\n");
        find_contacts(12.0, &pos[0], (Float *)NULL,
                A.begin(), A.end(), B.begin(), B.end(), Output());
        printf("  PERIODIC\n");
        find_contacts(12.0, &pos[0], mol->global_cell[0], 
                A.begin(), A.end(), B.begin(), B.end(), Output());
    }
    return 0;
}

