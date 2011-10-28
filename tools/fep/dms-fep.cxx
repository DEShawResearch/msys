#include "dms.hxx"
#include "alchemical.hxx"

#include <fstream>
#include <assert.h>
#include <stdio.h>

using namespace desres::msys;

static void read_map(std::istream& in, IdList& amap, IdList& bmap) {

    char buf[256];
    while (in.getline(buf, sizeof(buf))) {
        std::streamsize n=in.gcount();
        if (n==sizeof(buf)-1) {
            fprintf(stderr, "Line too long in atom map\n");
            exit(1);
        }
        if (n==1) continue;
        assert(n>1);
        if (buf[0]=='#') continue;
        int a,b;
        if (sscanf(buf, "%d %d", &a, &b)!=2) {
            fprintf(stderr, "Misformatted line in atom map: %s\n", buf);
            exit(1);
        }
        if (a==0 || b==0) {
            fprintf(stderr, "atom map contains an invalid zero: %s\n", buf);
            exit(1);
        }
        amap.push_back( a<0 ? BadId : a-1 );
        bmap.push_back( b<0 ? BadId : b-1 );
    }
}

int main(int argc, char *argv[]) {
    if (argc!=5) {
        fprintf(stderr, "Usage: %s A.dms B.dms atom.map C.dms\n", argv[0]);
        exit(1);
    }

    SystemPtr A = ImportDMS(argv[1]);
    SystemPtr B = ImportDMS(argv[2]);
    std::ifstream in(argv[3]);
    if (!in) {
        fprintf(stderr, "Unable to read atom map at path %s\n", argv[3]);
        exit(1);
    }
    IdList amap, bmap;
    read_map(in, amap, bmap);

    SystemPtr C = CreateAlchemical( A, amap, B, bmap );
    ExportDMS(C, argv[4]);
    return 0;
}

