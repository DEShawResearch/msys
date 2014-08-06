#include "dms.hxx"
#include "term_table.hxx"
#include <cstdio>
#include <boost/foreach.hpp>

using namespace desres::msys;

int main(int argc, char *argv[]) {
    if (argc!=3) {
        fprintf(stderr, "Usage: %s input.dms tablename\n",
                argv[0]);
        return 1;
    }
    SystemPtr sys = ImportDMS(argv[1]);
    TermTablePtr table = sys->table(argv[2]);
    for (TermTable::const_iterator i=table->begin(),e=table->end();i!=e; ++i) {
        printf("id %u atom %u %u param %u\n", 
                i->id(), i->atom(0), i->atom(1), i->param());
    }
    //TermTable const& T = *table;
    //BOOST_FOREACH(TermTable::term_t t, T) {
    //}
    return 0;
}

