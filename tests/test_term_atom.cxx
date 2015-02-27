#include "io.hxx"
#include <boost/foreach.hpp>
#include <stdio.h>

using namespace desres::msys;


int main(int argc, char *argv[]) {
    for (int i=1; i<argc; i++) {
        SystemPtr mol = Load(argv[i]);
        std::vector<std::string> names = mol->tableNames();

        const int ntries=5;
        std::vector<double> timings;

        printf("v2mdsim way: iterate terms, atoms, insert into set\n");
        timings.clear();
        timings.reserve(ntries);
        for (int k=0; k<ntries; k++) {
            double t=-now();
            std::set<Id> atoms;
            /* for each table */
            BOOST_FOREACH(std::string const& name, names) {
                TermTablePtr table = mol->table(name);
                BOOST_FOREACH(Id term, table->terms()) {
                    for (Id j=0; j<table->atomCount(); j++) {
                        atoms.insert(table->atom(term,j));
                    }
                }
            }
            t += now();
            printf("atoms: %lu\n", atoms.size());
            timings.push_back(t);
        }
        std::sort(timings.begin(), timings.end());
        printf("best: %lf worst %lf\n", timings.front(), timings.back());

        printf("v2mdsim, but append to vector, then sort_unique\n");
        timings.clear();
        timings.reserve(ntries);
        for (int k=0; k<ntries; k++) {
            double t=-now();
            IdList atoms;
            atoms.reserve(mol->maxAtomId());
            /* for each table */
            BOOST_FOREACH(std::string const& name, names) {
                TermTablePtr table = mol->table(name);
                BOOST_FOREACH(Id term, table->terms()) {
                    for (Id j=0; j<table->atomCount(); j++) {
                        atoms.push_back(table->atom(term,j));
                    }
                }
            }
            sort_unique(atoms);
            t += now();
            printf("atoms: %lu\n", atoms.size());
            timings.push_back(t);
        }
        std::sort(timings.begin(), timings.end());
        printf("best: %lf worst %lf\n", timings.front(), timings.back());

        printf("Use msys term iterator into vector\n");
        timings.clear();
        timings.reserve(ntries);
        for (int k=0; k<ntries; k++) {
            double t=-now();
            IdList atoms;
            atoms.reserve(mol->maxAtomId());
            /* for each table */
            BOOST_FOREACH(std::string const& name, names) {
                TermTablePtr table = mol->table(name);
                for (TermTable::const_iterator b=table->begin(),e=table->end(); b!=e; ++b) {
                    atoms.insert(atoms.end(), b->atoms(), b->atoms()+b->size());
                }
            }
            sort_unique(atoms);
            t += now();
            printf("atoms: %lu\n", atoms.size());
            timings.push_back(t);
        }
        std::sort(timings.begin(), timings.end());
        printf("best: %lf worst %lf\n", timings.front(), timings.back());
    }
    return 0;
}
