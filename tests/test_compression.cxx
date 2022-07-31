#include <msys/compression.hxx>
#include <fstream>

int main(int argc, char *argv[]) {
    static int BUFSIZE = 4096;
    char buf[BUFSIZE+1];
    for (int i=1; i<argc; i++) {
        std::ifstream in(argv[i]);
        auto iptr = desres::msys::maybe_compressed_istream(in);
        while (*iptr) {
            iptr->read(buf, BUFSIZE);
            buf[iptr->gcount()] = '\0';
            std::cout << buf;
        }
    }
    return 0;
}


