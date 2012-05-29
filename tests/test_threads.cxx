#include <msys/dms.hxx>
#include <boost/thread.hpp>
#include <boost/foreach.hpp>

using namespace desres::msys;

static void worker(std::string ifile) {
    printf("loading %s\n", ifile.c_str());
    ImportDMS(ifile);
    printf("finished %s\n", ifile.c_str());
}

int main(int argc, char *argv[]) {
    std::vector<boost::thread> v;
    for (int i=1; i<argc; i++) {
        v.push_back(boost::thread(worker, argv[i]));
    }
    BOOST_FOREACH(boost::thread& t, v) t.join();
    return 0;
}
