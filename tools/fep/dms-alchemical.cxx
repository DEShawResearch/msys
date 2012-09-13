#include "load.hxx"
#include "dms.hxx"
#include "alchemical.hxx"
#include <boost/algorithm/string.hpp>
#include <fstream>
#include <string.h>

using namespace desres::msys;

static bool avoid_alchemical_noop = true;

void parse_cmdline( int *pargc, char ***pargv ) {
    int i,j=0;
                           
    for (i=0; i<pargc[0]; i++) {
        pargv[0][j] = pargv[0][i];
        if ( !strcmp(pargv[0][j], "--keep-alchemical-noop")) {
            avoid_alchemical_noop = false; 
        }  else {
            j++;
        }
    }
    pargc[0]=j;
    pargv[0][j]=NULL;
}
static std::vector<IdPair> read_atommap(const char *path) {
    std::vector<IdPair> pairs;
    char buf[80];
    std::ifstream in(path);
    if (!in) {
        MSYS_FAIL("Unable to read atommap at '" << path << "'");
    }
    int lineno = 0;
    while (in.getline(buf,sizeof(buf))) {
        ++lineno;
        std::string line(buf);
        boost::trim(line);
        if (line.empty() || line[0]=='#') {
            continue;
        }
        int a,b;
        if (sscanf(line.c_str(), "%d %d", &a, &b)!=2) {
            MSYS_FAIL("Misformtted atom map at line " << lineno);
        }
        if (a==0 || b==0) {
            MSYS_FAIL("Atom map entry contains a zero at line " << lineno);
        }
        Id A = a<0 ? BadId : a-1;
        Id B = b<0 ? BadId : b-1;
        pairs.push_back(std::make_pair(A,B));
    }
    return pairs;
}

int main(int argc, char *argv[]) {
    parse_cmdline(&argc, &argv);
    if (argc!=5) {
        fprintf(stderr, "\nUsage: %s A.dms B.dms atom.map C.dms\n", argv[0]);
        exit(1);
    }

    SystemPtr A = Load(argv[1]);
    SystemPtr B = Load(argv[2]);
    std::vector<IdPair> pairs = read_atommap(argv[3]);
    SystemPtr C = MakeAlchemical(A,B,pairs,avoid_alchemical_noop);
    ExportDMS(C, argv[4], Provenance::fromArgs(argc,argv));
    return 0;
}

