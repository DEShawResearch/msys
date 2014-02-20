#include <stdio.h>
#include "cealign.hxx"

using namespace desres::msys;

static std::vector<float> load_positions(const char* path) {
    FILE* fp = fopen(path, "r");
    fseek(fp, 0, SEEK_END);
    long nbytes = ftell(fp);
    long nfloats = nbytes/4;
    assert(nfloats*4==nbytes);
    assert(((nfloats/3)*3)==nfloats);
    std::vector<float> pos(nfloats);
    fseek(fp, 0, SEEK_SET);
    fread(&pos[0], pos.size(), sizeof(float), fp);
    fclose(fp);

    printf("got %u positions from %s\n", (unsigned)pos.size()/3, path);
    for (unsigned i=0, n=pos.size()/3; i<n; i++) {
        break;
        printf("%3u %8.3f %8.3f %8.3f\n",
                i,
                pos[3*i+0],
                pos[3*i+1],
                pos[3*i+2]);
    }
    return pos;
}

int main(int argc, char *argv[]) {
    std::vector<float> apos = load_positions(argv[1]);
    std::vector<float> bpos = load_positions(argv[2]);
    unsigned nA = apos.size()/3;
    unsigned nB = bpos.size()/3;
    unsigned n = std::max(nA, nB);
    std::vector<unsigned> inds(n);
    for (unsigned i=0; i<n; i++) inds[i]=i;

    CEAlign<float> ce = CEAlign<float>::WithDefaults();
    float mat[9];
    double rms = ce.compute(nA, &inds[0], &apos[0],
                            nB, &inds[0], &bpos[0], mat);
    printf("rms: %le\n", rms);
    return 0;
}

