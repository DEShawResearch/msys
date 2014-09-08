#include "io.hxx"
#include "dms.hxx"
#include "alchemical.hxx"
#include <boost/algorithm/string.hpp>
#include <fstream>
#include <string.h>

using namespace desres::msys;

static bool avoid_alchemical_noop = true;
static bool keep_none = false;

static void usage(FILE *fd) {
    fprintf(fd, "Usage:\ndms-alchemical A.dms B.dms atom.map C.dms\n");
    fprintf(fd,

"Create an alchemical system from A and B states and a map between them.\n"
"\n"
"The atom.map file should consist of lines with two 1-based indices,\n"
"the first referring to atoms in the A state and the second to atoms in\n"
"the B state.  Either the A or B index may be negative, indicating that\n"
"the corresponding atom has no analog in the other state.  The order of\n"
"the lines in the file is insignificant.\n"
"\n"
"It is not necessary that the atom map reference every atom in A and B\n"
"states; however, any term in a given state must be either completely\n"
"mapped or completely unmapped.  In practice this usually means that\n"
"the atom map should contain complete sets of connected atoms.\n"
"\n"
"Options:\n"
"   -h, --help       show this help message and exit\n"
"   --keep-alchemical-noop\n"
"                    Generate all alchemical terms described by the atom\n"
"                    map, even those whose A and B states are identical.\n"
"                    This option is present only for comparison with\n"
"                    previous versions of dms-alchemical.\n"
"   --keep-none\n"
"                    Do not keep any interactions between real atoms and\n"
"                    dummy atoms.\n"
    );
}

void parse_cmdline( int *pargc, char ***pargv ) {
    int i,j=0;
                           
    for (i=0; i<pargc[0]; i++) {
        pargv[0][j] = pargv[0][i];
        if ( !strcmp(pargv[0][j], "-h") ||
             !strcmp(pargv[0][j], "--help")) {
            usage(stdout);
            exit(0);

        } else if ( !strcmp(pargv[0][j], "--keep-alchemical-noop")) {
            avoid_alchemical_noop = false; 
        } else if ( !strcmp(pargv[0][j], "--keep-none")) {
            keep_none = true;
        } else {
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
        usage(stderr);
        exit(1);
    }

    SystemPtr A = Load(argv[1]);
    SystemPtr B = Load(argv[2]);
    std::vector<IdPair> pairs = read_atommap(argv[3]);
    SystemPtr C = MakeAlchemical(A,B,pairs,avoid_alchemical_noop,keep_none);
    ExportDMS(C, argv[4], Provenance::fromArgs(argc,argv));
    return 0;
}

