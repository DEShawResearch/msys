#include "atomsel.hxx"
#include "io.hxx"
#include "pfx/pfx.hxx"
#include "spatial_hash.hxx"

#include <stdlib.h>

using namespace desres::msys;

int main(int argc, char *argv[]) {

    /* 
     * generate N random points bounded by [0,1).
     * Find the number of points pbwithin 0.1 of the first point.
     * Rotate box and cell by a random rotation.
     * Check that number of points remains the same.
     */

    static const int N = 100;
    static const float radius = 0.7;
    IdList wsel(N), psel(1,0);
    srand48(1973);
    double cell[9] = {2,0,0, 0,3,0, 0,0,4};
    float pos[3*N];
    for (int i=0; i<N; i++) {
        pos[3*i  ] = drand48()*cell[0];
        pos[3*i+1] = drand48()*cell[4];
        pos[3*i+2] = drand48()*cell[8];
        wsel[i] = i;
    }
    IdList ids;

    ids = FindWithin(wsel, pos, psel, pos, radius, cell);
    printf("Found %lu\n", ids.size());
    for (unsigned i=0, n=ids.size(); i<n; i++) {
        printf("%2u ", ids[i]);
    }
    printf("\n");

    /* apply rotation */
    float R[9];
    memset(R,0,sizeof(R));
    float iroot2 = 1/sqrt(2);
    R[0]= iroot2;
    R[1]=-iroot2;
    R[3]= iroot2;
    R[4]= iroot2;
    R[8]=1;
    pfx::apply_rotation(N,pos,R);
    {
        float box[9];
        std::copy(cell,cell+9,box);
        pfx::apply_rotation(3,box,R);
        std::copy(box,box+9,cell);
    }

    ids = FindWithin(wsel, pos, psel, pos, radius, cell);
    printf("Found %lu\n", ids.size());
    for (unsigned i=0, n=ids.size(); i<n; i++) {
        printf("%2u ", ids[i]);
    }
    printf("\n");


    return 0;
}
