#include <msys/dms.hxx>
#include <boost/thread.hpp>
#include <boost/foreach.hpp>
#include <cstdio>

using namespace desres::msys;

static void worker(std::string ifile) {
    printf("loading %s\n", ifile.c_str());
    ImportDMS(ifile);
    printf("finished %s\n", ifile.c_str());
}

typedef boost::shared_ptr<boost::thread> ThreadPtr;

int main(int argc, char *argv[]) {
    std::vector<ThreadPtr> v;
    for (int i=1; i<argc; i++) {
        v.push_back(ThreadPtr(new boost::thread(worker, argv[i])));
    }
    BOOST_FOREACH(ThreadPtr t, v) t->join();
    return 0;
}
