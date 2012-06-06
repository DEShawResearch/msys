#include "system.hxx"
#include "dms.hxx"

#include <sys/time.h>

using namespace desres::msys;

static double time_of_day() {
  struct timeval tm;
  struct timezone tz;

  gettimeofday(&tm, &tz);
  return((double)(tm.tv_sec) + (double)(tm.tv_usec)/1000000.0);
}

static SystemPtr make_system(int n) {
    SystemPtr m = System::create();
    Id chn = m->addChain();
    for (int i=0; i<n; i++) {
        Id res = m->addResidue(chn);
        m->residue(res).resid = i;
        m->residue(res).name = "TIP3";

        Id atm = m->addAtom(res);
        m->atom(atm).atomic_number = 8;
        m->atom(atm).x = drand48();
        m->atom(atm).y = drand48();
        m->atom(atm).z = drand48();
        m->atom(atm).vx = drand48();
        m->atom(atm).vy = drand48();
        m->atom(atm).vz = drand48();
        m->atom(atm).mass = 16;
        m->atom(atm).name = "OH";

        atm = m->addAtom(res);
        m->atom(atm).atomic_number = 1;
        m->atom(atm).x = drand48();
        m->atom(atm).y = drand48();
        m->atom(atm).z = drand48();
        m->atom(atm).vx = drand48();
        m->atom(atm).vy = drand48();
        m->atom(atm).vz = drand48();
        m->atom(atm).mass = 1;
        m->atom(atm).name = "H1";

        atm = m->addAtom(res);
        m->atom(atm).atomic_number = 1;
        m->atom(atm).x = drand48();
        m->atom(atm).y = drand48();
        m->atom(atm).z = drand48();
        m->atom(atm).vx = drand48();
        m->atom(atm).vy = drand48();
        m->atom(atm).vz = drand48();
        m->atom(atm).mass = 1;
        m->atom(atm).name = "H2";
    }

    return m;
}


int main(int argc, char *argv[]) {
    printf("%8s %12s %12s %12s %12s\n", "ATOMS", "MAKE", "SAVE", "LOAD",
            "LOAD/MAKE");
    for (int i=1; i<argc; i++) {
        int n = atoi(argv[i]);
        double t1=time_of_day();
        SystemPtr m = make_system(n);
        double t2 = time_of_day();
        ExportDMS(m, "/tmp/msys_bench.dms", Provenance::fromArgs(argc,argv));
        m.reset();
        double t3 = time_of_day();
        m = ImportDMS("/tmp/msys_bench.dms");
        double t4 = time_of_day();
        printf("%8d %12f %12f %12f %12f\n", n, t2-t1, t3-t2, t4-t3,
                (t4-t3)/(t2-t1));
    }
    return 0;
}
