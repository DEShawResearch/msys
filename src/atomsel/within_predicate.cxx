/* @COPYRIGHT@ */

#include "within_predicate.hxx"
#include "msys_keyword.hxx"
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <algorithm>
#include <limits>


#define STACK_SIZE 4

using namespace desres::msys::atomsel;

typedef desres::msys::atom_t atom;
using desres::msys::SystemPtr;
using desres::msys::atom_t;
using desres::msys::now;
using desres::msys::Id;
using desres::msys::IdSet;
using desres::msys::IdList;

struct point_t {
    double x,y,z;
    point_t() {}
    point_t(double _x, double _y, double _z) : x(_x), y(_y), z(_z) {}
};

struct posarray : boost::noncopyable {
    std::vector<point_t> pos;
    double xmin, ymin, zmin;
    double xmax, ymax, zmax;
    Id size;

    posarray(SystemPtr mol, Selection const& S)
    : xmin(), ymin(), zmin(), xmax(), ymax(), zmax(), size(S.count()) {
        pos.resize(size);
        for (Id i=0, j=0, n=S.size(); i<n; i++) {
            if (S[i]) {
                double x = mol->atomFAST(i).x;
                double y = mol->atomFAST(i).y;
                double z = mol->atomFAST(i).z;
                if (j==0) {
                    xmin = xmax = x;
                    ymin = ymax = y;
                    zmin = zmax = z;
                } else {
                    xmin = std::min(xmin, x);
                    ymin = std::min(ymin, y);
                    zmin = std::min(zmin, z);
                    xmax = std::max(xmax, x);
                    ymax = std::max(ymax, y);
                    zmax = std::max(zmax, z);
                }
                pos[j].x = x;
                pos[j].y = y;
                pos[j].z = z;
                ++j;
            }
        }
    }
    double const& x(Id i) const { return pos[i].x; }
    double const& y(Id i) const { return pos[i].y; }
    double const& z(Id i) const { return pos[i].z; }
};


struct voxel_t {
    /* For voxels on the edge of the grid, compute and store the
     * indices of neighbor voxels.  For other voxels, the neighbors
     * are computed using a cached set of offsets. */
  int* nbrs;  /* NULL for non-edge voxel, terminated by -1 */
  Id stack[STACK_SIZE];
  Id * points;
  int num;
  int max;

  voxel_t() : nbrs(), points(stack), num(), max(STACK_SIZE) {}
  ~voxel_t() {
      if (points!=stack) free(points);
      if (nbrs) free(nbrs);
  }

  void add(Id p) {
      if (num==max) {
          max *= 1.3;
          if (points==stack) {
              points = (Id*)malloc(max*sizeof(Id));
              memcpy(points, stack, num*sizeof(Id));
          } else {
              points = (Id*)realloc(points, max*sizeof(Id));
          }
      }
      points[num++] = p;
  }
};

static
void find_voxel_full_shell_neighbors(voxel_t *mesh, int nx, int ny, int nz) {
  for (int zi=0; zi<nz; zi++) {
    for (int yi=0; yi<ny; yi++) {
      for (int xi=0; xi<nx; xi++) {
        if (zi!=0 && zi!=nz-1 &&
            yi!=0 && yi!=ny-1 &&
            xi!=0 && xi!=nx-1) {
            continue;
        }
        int n=0;
        int* nbrs = (int*)malloc(28*sizeof(int));
        int self = xi + nx*(yi + ny*zi);
        // it's a big win to always search the self voxel first!
        if (mesh[self].num) nbrs[++n]=self;
        for (int ti=zi-1; ti<=zi+1; ti++) {
          if (ti<0 || ti>=nz) continue;
          for (int si=yi-1; si<=yi+1; si++) {
            if (si<0 || si>=ny) continue;
            for (int ri=xi-1; ri<=xi+1; ri++) {
              if (ri<0 || ri>=nx) continue;
              int index = ri + nx*(si + ny*ti);
              if (index!=self && mesh[index].num) nbrs[++n] = index;
            }
          }
        }
        nbrs[0] = n;
        mesh[self].nbrs = nbrs;
      }
    }
  }
}

static
void compute_central_voxel_neighbors(int *nbrs, int nx, int ny, int nz) {
    *nbrs++ = 0;
    for (int k=-1; k<=1; k++) {
        for (int j=-1; j<=1; j++) {
            for (int i=-1; i<=1; i++) {
                if (!(k==0 && j==0 && i==0)) {
                    *nbrs++ = i + nx*(j + ny*k);
                }
            }
        }
    }
}

namespace {
  using desres::msys::SystemPtr;
  class WithinPredicate : public Predicate {
    SystemPtr sys;
    const double rad;
    PredicatePtr sub;
    const bool exclude;
    const bool periodic;

    public:
    WithinPredicate( SystemPtr e, double r, bool excl, bool per, PredicatePtr s )
      : sys(e), rad(r), sub(s), exclude(excl), periodic(per) {}

    void eval( Selection& s );
    void dump( std::ostream& str ) const {
      if (periodic) str << "pb";
      if (exclude) str << "ex";
      str << "within " << rad << " of [";
      sub->dump(str);
      str << "]";
    }
  };
  class WithinBondsPredicate : public Predicate {
    SystemPtr sys;
    const int N;
    PredicatePtr sub;

    public:
    WithinBondsPredicate( SystemPtr e, int n, PredicatePtr s )
      : sys(e), N(n), sub(s) {}

    void eval( Selection& s );
    void dump( std::ostream& str ) const {
      str << "withinbonds " << N << " of [";
      sub->dump(str);
      str << "]";
    }
  };
}

static void find_within( const point_t* wat,
                         posarray const& pro,
                         Selection& S,
                         Selection const& subsel,
                         double rad,
                         bool exclude ) {

  double t1=1000*now();

  /* If this is an "exwithin" selection, AND S with the converse of the
   * subselection */
  if (exclude) {
    S.subtract(subsel);
  }

  if (rad <= 0) {
    S.intersect(subsel);
    return;
  }
  const double ir = 1.0/rad;
  const double r2 = rad*rad;

  /* extend bounding box by selection radius */
  double xmin = pro.xmin - rad;
  double ymin = pro.ymin - rad;
  double zmin = pro.zmin - rad;
  double xsize = rad + pro.xmax - xmin;
  double ysize = rad + pro.ymax - ymin;
  double zsize = rad + pro.zmax - zmin;

  /* create and initialize voxel mesh */
  int nx = (int)(xsize/rad)+1;
  int ny = (int)(ysize/rad)+1;
  int nz = (int)(zsize/rad)+1;
  int nvoxel = nx*ny*nz;

  double t2=1000*now();

  /* map pro atoms to voxels */
  voxel_t * mesh = new voxel_t[nvoxel];
  for (Id i=0, n=pro.size; i<n; i++) {
    int xi = (pro.x(i)-xmin)*ir;
    int yi = (pro.y(i)-ymin)*ir;
    int zi = (pro.z(i)-zmin)*ir;
    int index = xi + nx*(yi + ny*zi);
    mesh[index].add(i);
  }

  double t3=1000*now();

  find_voxel_full_shell_neighbors(mesh, nx, ny, nz);
  int central_nbrs[27];
  compute_central_voxel_neighbors(central_nbrs, nx, ny, nz);

  double t4=1000*now();

  /* loop over wat atoms */
  for (Id i=0; i<S.size(); i++) {
    if (!S[i]) continue;
    if (subsel[i]) continue;
    point_t const& p = wat[i];
    const double x = p.x;
    const double y = p.y;
    const double z = p.z;
    int xi = (x-xmin)*ir;
    int yi = (y-ymin)*ir;
    int zi = (z-zmin)*ir;
    if (xi<0 || xi>=nx ||
        yi<0 || yi>=ny ||
        zi<0 || zi>=nz) {
      S[i]=0;
      continue;
    }
    int index = xi + nx*(yi + ny*zi);
    const voxel_t * v = mesh+index;
    const int * nbrs;
    int n_nbrs;
    int self;
    if (v->nbrs) {
        /* edge voxel, use existing neighbors */
        nbrs = v->nbrs;
        n_nbrs = *nbrs++;
        self = 0;   /* absolute indices */
    } else {
        nbrs = central_nbrs;
        n_nbrs = 27;
        self = index;
    }
    bool on=0;
    for (int j=0; j<n_nbrs; j++) {
      const voxel_t* nbr = mesh + nbrs[j]+self;
      int natoms = nbr->num;
      for (int k=0; k<natoms; k++) {
        const Id pk = nbr->points[k];
        double dx=x-pro.x(pk);
        double dy=y-pro.y(pk);
        double dz=z-pro.z(pk);
        if (dx*dx + dy*dy + dz*dz <=r2) {
          on=1;
          break;
        }
      }
      if (on) break;
    }
    if (!on) S[i]=0;
  }
  double t5=1000*now();

  delete [] mesh;

  double t6=1000*now();
  if (0) printf("within %8.3f: tot %8.3f : %f %f %f %f %f\n", rad,t6-t1,
          t2-t1,t3-t2,t4-t3,t5-t4,t6-t5);
}

void WithinPredicate::eval( Selection& S ) {
    Selection subsel = full_selection(sys);
    sub->eval(subsel);

    if (rad<=0) {
        if (exclude) S.subtract(subsel);
        S.intersect(subsel);
        return;
    }

    /* "protein" coordinates */
    posarray pro(sys, subsel);
    if (pro.size==0) {
        S.clear();
        return;
    }

    /* "water" coordinates */
    std::vector<point_t> wat(S.size());
    for (unsigned i=0; i<S.size(); i++) {
        if (S[i]) {
            wat[i].x = sys->atomFAST(i).x;
            wat[i].y = sys->atomFAST(i).y;
            wat[i].z = sys->atomFAST(i).z;
        }
    }

    if (periodic) {
        /* replicate the water within a bounding box around protein */
        double xmin = pro.xmin - rad;
        double ymin = pro.ymin - rad;
        double zmin = pro.zmin - rad;
        double xmax = pro.xmax + rad;
        double ymax = pro.ymax + rad;
        double zmax = pro.zmax + rad;

        IdList repids;  /* ids of replicated atoms */
        const double* A = sys->global_cell[0];
        const double* B = sys->global_cell[1];
        const double* C = sys->global_cell[2];
        if (A[1] || A[2] || 
            B[0] || B[2] ||
            C[0] || C[1]) {
            MSYS_FAIL("pbwithin does not support triclinic global cell");
        }
        double ga = A[0];
        double gb = B[1];
        double gc = C[2];
        for (Id id=0; id<S.size(); id++) {
            if (!S[id]) continue;
            /* Need a copy since we'll append to wat below */
            point_t p = wat[id];    
            for (int i=-1; i<=1; i++) {
                double x = p.x + ga*i;
                if (x<xmin || x >= xmax) continue;
                for (int j=-1; j<=1; j++) {
                    double y = p.y + gb*j;
                    if (y<ymin || y >= ymax) continue;
                    for (int k=-1; k<=1; k++) {
                        double z = p.z + gc*k;
                        if (z<zmin || z >= zmax) continue;
                        if (i==0 && j==0 && k==0) continue;
                        repids.push_back(id);
                        wat.push_back(point_t(x,y,z));
                    }
                }
            }
        }

        Selection Srep(wat.size());
        Selection Ssub(wat.size());
        /* copy original flags */
        for (unsigned i=0; i<S.size(); i++) Srep[i]=S[i];
        /* turn on all replica flags - we know they're selected */
        for (unsigned i=0; i<repids.size(); i++) Srep[i+S.size()]=1;
        /* copy subsel into expanded set */
        for (unsigned i=0; i<subsel.size(); i++) Ssub[i]=subsel[i];

        find_within( &wat[0], pro, Srep, Ssub, rad, false );

        /* copy Srep back into S */
        for (unsigned i=0; i<S.size(); i++) S[i]=Srep[i];

        /* OR the replica flags into the home set */
        for (Id i=0; i<repids.size(); i++) {
            if (Srep[S.size()+i]) {
                S[repids[i]]=1;
            }
        }

    } else { 
        find_within( &wat[0], pro, S, subsel, rad, exclude );
    }
}

void WithinBondsPredicate::eval( Selection& S ) {
  Selection subsel = full_selection(sys);
  sub->eval(subsel);
  IdSet atms;
  /* hash the atoms in the subsel */
  for (Id i=0; i<subsel.size(); i++) {
    if (subsel[i]) atms.insert(i);
  }
  /* expand subsel by N bonds */
  for (int i=0; i<N; i++) {
    IdList tmp;
    for (IdSet::const_iterator iter=atms.begin(); iter!=atms.end(); ++iter) {
        IdList bonds = sys->bondsForAtom(*iter);
        for (Id j=0; j<bonds.size(); j++) {
            tmp.push_back(sys->bond(bonds[j]).other(*iter));
        }
    }
    atms.insert(tmp.begin(), tmp.end());
  }
  for (Id i=0; i<subsel.size(); i++) {
    subsel[i] = atms.count(i);
  }
  S.intersect(subsel);
}

namespace desres { namespace msys { namespace atomsel {

  PredicatePtr within_predicate( SystemPtr sys, double r, PredicatePtr s ) {
    return PredicatePtr(new WithinPredicate(sys,r,false,false,s));
  }

  PredicatePtr exwithin_predicate( SystemPtr sys, double r, PredicatePtr s ) {
    return PredicatePtr(new WithinPredicate(sys,r,true,false,s));
  }

  PredicatePtr pbwithin_predicate( SystemPtr sys, double r, PredicatePtr s ) {
    return PredicatePtr(new WithinPredicate(sys,r,false,true,s));
  }

  PredicatePtr withinbonds_predicate( SystemPtr sys, int n, PredicatePtr s ) {
    return PredicatePtr(new WithinBondsPredicate(sys,n,s));
  }

}}}

namespace {
    class KNearestPredicate : public Predicate {
        SystemPtr _sys;
        const unsigned _N;
        PredicatePtr _sub;

    public:
        KNearestPredicate(SystemPtr sys, unsigned k, PredicatePtr sub)
            : _sys(sys), _N(k), _sub(sub) {}

        void eval(Selection& s);
        void dump( std::ostream& str ) const {
            str << "nearest " << _N << " to [";
            _sub->dump(str);
            str << "]";
        }
    };
}

void KNearestPredicate::eval( Selection& S ) {

    /* evaluate subselection */
    Selection subsel = full_selection(_sys);
    _sub->eval(subsel);

    /* "protein" coordinates */
    posarray pro(_sys, subsel);
    if (pro.size==0) {
        S.clear();
        return;
    }

    /* we never include subselection */
    S.subtract(subsel);

    /* more than k atoms available? */
    double rmin = 0;
    Selection smin(S);
    Id nmin = smin.count();
    if (nmin <= _N) return;

    /* "water" coordinates */
    std::vector<point_t> wat(S.size());
    for (unsigned i=0; i<S.size(); i++) {
        if (S[i]) {
            wat[i].x = _sys->atomFAST(i).x;
            wat[i].y = _sys->atomFAST(i).y;
            wat[i].z = _sys->atomFAST(i).z;
        }
    }

    double rmax=2.5;
    Selection smax(0);
    Id nmax=nmin;

    /* increase rmax until Ny is at least _N */
    for (;;) {
        smax = S;
        find_within( &wat[0], pro, smax, subsel, rmax, false);
        nmax = smax.count();
        //printf("rmin %f nmin %u rmax %f nmax %u\n", rmin, nmin, rmax, nmax);
        if (nmax >= _N) break;
        rmin = rmax;
        smin = smax;
        nmin = nmax;
        rmax *= 1.5;
    }

    /* Do a couple rounds of bisection search to narrow it down */
    for (int nb=0; nb<6; nb++) {
        Selection sm(S);
        double rm = 0.5*(rmin+rmax);
        find_within( &wat[0], pro, sm, subsel, rm, false);
        Id nm = sm.count();
        //printf("rm %f nm %u\n", rm, nm);
        if (nm>=_N) {
            smax = sm;
            rmax = rm;
            nmax = nm;
        } else {
            smin = sm;
            rmin = rm;
            nmin = nm;
        }
    }
    //printf("min: rad %f n %u\n", rmin, nmin);
    //printf("max: rad %f n %u\n", rmax, nmax);

    std::vector<std::pair<double,Id> > pts;
    /* for each water in smax but not in smin */
    for (Id i=0, n=S.size(); i<n; i++) {
        if (smin[i] || !smax[i]) continue;
        point_t const& w = wat[i];
        double r2 = std::numeric_limits<double>::max();
        /* find min dist to protein */
        for (Id j=0, m=pro.size; j<m; j++) {
            double dx = w.x-pro.x(j);
            double dy = w.y-pro.y(j);
            double dz = w.z-pro.z(j);
            double d2 = dx*dx + dy*dy + dz*dz;
            r2 = std::min(r2, d2);
        }
        pts.push_back(std::make_pair(r2, i));
    }
    std::partial_sort(pts.begin(), pts.begin()+(_N-nmin), pts.end());
    S = smin;
    for (Id i=0, n=_N-nmin; i<n; i++) {
        S[pts[i].second] = 1;
    }
}

namespace desres { namespace msys { namespace atomsel {

    PredicatePtr k_nearest_predicate(SystemPtr sys, unsigned k, PredicatePtr S) {
        return PredicatePtr(new KNearestPredicate(sys,k,S));
    }

}}}


