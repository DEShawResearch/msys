#include "system.hxx"
#include "dms.hxx"
#include "schema.hxx"

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
    TermTablePtr stretch = AddTable(m, "stretch_harm");
    TermTablePtr angle = AddTable(m, "angle_harm");
    TermTablePtr excl = AddTable(m, "exclusion");
    TermTablePtr cons = AddTable(m, "constraint_hoh");

    stretch->params()->addParam();
    angle->params()->addParam();
    cons->params()->addParam();

    IdList pair(2), triple(3);

    for (int i=0; i<n; i++) {
        Id res = m->addResidue(chn);
        m->residue(res).resid = i;
        m->residue(res).name = "TIP3";

        Id atm;

        Id OH = atm = m->addAtom(res);
        m->atom(atm).atomic_number = 8;
        m->atom(atm).x = drand48();
        m->atom(atm).y = drand48();
        m->atom(atm).z = drand48();
        m->atom(atm).vx = drand48();
        m->atom(atm).vy = drand48();
        m->atom(atm).vz = drand48();
        m->atom(atm).mass = 16;
        m->atom(atm).name = "OH";

        Id H1 = atm = m->addAtom(res);
        m->atom(atm).atomic_number = 1;
        m->atom(atm).x = drand48();
        m->atom(atm).y = drand48();
        m->atom(atm).z = drand48();
        m->atom(atm).vx = drand48();
        m->atom(atm).vy = drand48();
        m->atom(atm).vz = drand48();
        m->atom(atm).mass = 1;
        m->atom(atm).name = "H1";

        Id H2 = atm = m->addAtom(res);
        m->atom(atm).atomic_number = 1;
        m->atom(atm).x = drand48();
        m->atom(atm).y = drand48();
        m->atom(atm).z = drand48();
        m->atom(atm).vx = drand48();
        m->atom(atm).vy = drand48();
        m->atom(atm).vz = drand48();
        m->atom(atm).mass = 1;
        m->atom(atm).name = "H2";

        m->addBond(OH,H1);
        m->addBond(OH,H2);

        pair[0] = OH;
        pair[1] = H1;
        stretch->addTerm(pair,0);
        excl->addTerm(pair,BadId);
        pair[2] = H2;
        stretch->addTerm(pair,0);
        excl->addTerm(pair,BadId);

        triple[0] = OH;
        triple[1] = H1;
        triple[2] = H2;
        angle->addTerm(triple,0);
        cons->addTerm(triple,0);
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
