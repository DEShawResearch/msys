/* @COPYRIGHT@ */

#include "within_predicate.hxx"
#include "msys_keyword.hxx"
#include <math.h>
#include <string.h>
#include <stdlib.h>

#define STACK_SIZE 16

using namespace desres::msys::atomsel;

typedef desres::msys::atom_t atom;
using desres::msys::SystemPtr;
using desres::msys::Id;
using desres::msys::IdSet;
using desres::msys::IdList;

struct point_t {
    double x,y,z;
};

struct voxel_t {
  int nbrs[27];
  int n_nbrs;
  point_t stack[STACK_SIZE];
  point_t * points;
  int num;
  int max;

  voxel_t() : n_nbrs(0), points(stack), num(0), max(STACK_SIZE) {}
  ~voxel_t() {
      if (points!=stack) free(points);
  }

  void add(point_t const& p) {
      if (num==max) {
          max *= 1.3;
          if (points==stack) {
              points = (point_t *)malloc(max*sizeof(point_t));
              memcpy(points, stack, num*sizeof(point_t));
          } else {
              points = (point_t *)realloc(points, max*sizeof(point_t));
          }
      }
      points[num++] = p;
  }
};


/* find the bounding box of selected atoms. Return 0 if no atoms 
 * selected. */
static
int find_bbox(const Selection& S, const point_t* points, 
              double *min, double *max) {
  double xmin, ymin, zmin;
  double xmax, ymax, zmax;

  // find first selected atom
  int i,n=S.size();
  for (i=0; i<n; i++) if (S[i]) break;
  if (i==n) return 0;

  // initialize min/max to first atom position
  xmin=xmax=points[i].x;
  ymin=ymax=points[i].y;
  zmin=zmax=points[i].z;

  // extend minmax
  for (;i<n; i++) {
    if (!S[i]) continue;
    point_t const& p = points[i];
    double x=p.x;
    double y=p.y;
    double z=p.z;

    if (x<xmin) xmin=x;
    if (x>xmax) xmax=x;

    if (y<ymin) ymin=y;
    if (y>ymax) ymax=y;

    if (z<zmin) zmin=z;
    if (z>zmax) zmax=z;
  }
  // store to supplied buffers
  min[0]=xmin;
  min[1]=ymin;
  min[2]=zmin;
  max[0]=xmax;
  max[1]=ymax;
  max[2]=zmax;
  return 1;
}

static
void find_voxel_full_shell_neighbors(voxel_t *mesh, int nx, int ny, int nz) {
  int i, zi, yi, xi, ti, si, ri;
  for (zi=0; zi<nz; zi++) {
    for (yi=0; yi<ny; yi++) {
      for (xi=0; xi<nx; xi++) {
        int n=0;
        int nbrs[27];
        int self = xi + nx*(yi + ny*zi);
        // it's a big win to always search the self voxel first!
        nbrs[n++]=self;
        for (ti=zi-1; ti<=zi+1; ti++) {
          if (ti<0 || ti>=nz) continue;
          for (si=yi-1; si<=yi+1; si++) {
            if (si<0 || si>=ny) continue;
            for (ri=xi-1; ri<=xi+1; ri++) {
              if (ri<0 || ri>=nx) continue;
              int index = ri + nx*(si + ny*ti);
              if (index!=self) nbrs[n++] = index;
            }
          }
        }
        mesh[self].n_nbrs=n;
        for (i=0; i<n; i++) mesh[self].nbrs[i] = nbrs[i];
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

static void find_within( const point_t* points,
                         Selection& S,
                         Selection const& subsel,
                         double rad,
                         bool exclude ) {

  double min[3], max[3];
  double xmin, ymin, zmin;
  double xsize, ysize, zsize;
  const double r2 = rad*rad;
  double ir;
  voxel_t * mesh;

  int nx, ny, nz, nvoxel;
  int i;

  /* find bounding box of subselection */
  if (!find_bbox(subsel, points, min, max)) {
    /* no atoms in subsel */
    S.clear();
    return;
  }

  /* If this is an "exwithin" selection, AND S with the converse of the
   * subselection */
  if (exclude) {
    S.subtract(subsel);
  }

  if (rad <= 0) {
    S.intersect(subsel);
    return;
  }
  ir = 1.0/rad;

  /* extend bounding box by selection radius */
  for (i=0; i<3; i++) {
    min[i] -= rad;
    max[i] += rad;
  }
  xmin=min[0];
  ymin=min[1];
  zmin=min[2];
  xsize=max[0]-xmin;
  ysize=max[1]-ymin;
  zsize=max[2]-zmin;

  /* create and initialize voxel mesh */
  nx = (int)(xsize/rad)+3;
  ny = (int)(ysize/rad)+3;
  nz = (int)(zsize/rad)+3;
  nvoxel = nx*ny*nz;

  mesh = new voxel_t[nvoxel];

  find_voxel_full_shell_neighbors(mesh, nx, ny, nz);

  /* map subsel atoms to voxels */
  for (Id i=0; i<subsel.size(); i++) {
    if (!subsel[i]) continue;

    point_t const& p = points[i];
    int xi = (p.x-xmin)*ir;
    int yi = (p.y-ymin)*ir;
    int zi = (p.z-zmin)*ir;
    int index = xi + nx*(yi + ny*zi);
    mesh[index].add(p);
  }

  /* loop over atoms in left selection */
  for (Id i=0; i<S.size(); i++) {
    int xi,yi,zi;
    int j, index;
    int n_nbrs;
    const int * nbrs;
    int on;
    voxel_t * v;
    if (!S[i]) continue;
    if (subsel[i]) continue;
    point_t const& p = points[i];
    xi = (p.x-xmin)*ir;
    yi = (p.y-ymin)*ir;
    zi = (p.z-zmin)*ir;
    if (xi<0 || xi>=nx ||
        yi<0 || yi>=ny ||
        zi<0 || zi>=nz) {
      S[i]=0;
      continue;
    }
    index = xi + nx*(yi + ny*zi);
    v = mesh+index;
    n_nbrs = v->n_nbrs;
    nbrs = v->nbrs;
    on=0;
    for (j=0; j<n_nbrs; j++) {
      const voxel_t* nbr = mesh + nbrs[j];
      int k, natoms = nbr->num;
      for (k=0; k<natoms; k++) {
        point_t const& pk = nbr->points[k];
        double dx=pk.x - p.x;
        double dy=pk.y - p.y;
        double dz=pk.z - p.z;
        if (dx*dx + dy*dy + dz*dz <=r2) {
          on=1;
          break;
        }
      }
      if (on) break;
    }
    if (!on) S[i]=0;
  }
  delete [] mesh;
}

void WithinPredicate::eval( Selection& S ) {
    Selection subsel = full_selection(sys);
    sub->eval(subsel);

    std::vector<point_t> points(S.size());
    for (unsigned i=0; i<S.size(); i++) {
        if (S[i] || subsel[i]) {
            points[i].x = sys->atom(i).x;
            points[i].y = sys->atom(i).y;
            points[i].z = sys->atom(i).z;
        }
    }
    if (!points.size()) {
        S.clear();
        return;
    }

    if (periodic) {
        /* replicate the atoms in S within a bounding box around subsel */
        double min[3], max[3];
        if (!find_bbox(subsel, &points[0], min, max)) {
            S.clear();
            return;
        }
        for (int i=0; i<3; i++) {
            min[i] -= rad;
            max[i] += rad;
        }
        IdList repids;  /* ids of replicated atoms */
        desres::msys::Vec3 const& A = sys->global_cell.A;
        desres::msys::Vec3 const& B = sys->global_cell.B;
        desres::msys::Vec3 const& C = sys->global_cell.C;
        for (int i=-1; i<=1; i++) {
            double ax = A[0]*i;
            double ay = A[1]*i;
            double az = A[2]*i;
            for (int j=-1; j<=1; j++) {
                double bx = B[0]*j;
                double by = B[1]*j;
                double bz = B[2]*j;
                for (int k=-1; k<=1; k++) {
                    if (i==0 && j==0 && k==0) continue;
                    double vx = C[0]*k + ax + bx;
                    double vy = C[1]*k + ay + by;
                    double vz = C[2]*k + az + bz;

                    for (Id id=0; id<S.size(); id++) {
                        if (!S[id]) continue;
                        point_t p = points[id];
                        p.x += vx;
                        p.y += vy;
                        p.z += vz;
                        if (p.x > min[0] && p.x <= max[0] &&
                            p.y > min[1] && p.y <= max[1] &&
                            p.z > min[2] && p.z <= max[2]) {
                            points.push_back(p);
                            repids.push_back(id);
                        }
                    }
                }
            }
        }

        Selection Srep(points.size());
        Selection Ssub(points.size());
        /* copy original flags */
        for (unsigned i=0; i<S.size(); i++) Srep[i]=S[i];
        /* turn on all replica flags - we know they're selected */
        for (unsigned i=0; i<repids.size(); i++) Srep[i+S.size()]=1;
        /* copy subsel into expanded set */
        for (unsigned i=0; i<subsel.size(); i++) Ssub[i]=subsel[i];

        find_within( &points[0], Srep, Ssub, rad, false );

        /* copy Srep back into S */
        for (unsigned i=0; i<S.size(); i++) S[i]=Srep[i];

        /* OR the replica flags into the home set */
        for (Id i=0; i<repids.size(); i++) {
            if (Srep[S.size()+i]) {
                S[repids[i]]=1;
            }
        }

    } else { 
        find_within( &points[0], S, subsel, rad, exclude );
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

