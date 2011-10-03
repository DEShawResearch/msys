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

typedef struct voxel_t {
  int nbrs[27];
  int n_nbrs;
  atom * stack_atoms[STACK_SIZE];
  atom ** atoms;
  int num_atoms;
  int max_atoms;
} voxel;

static inline
void voxel_init( voxel * v ) {
  v->n_nbrs = 0;
  v->atoms = v->stack_atoms;
  v->num_atoms = 0;
  v->max_atoms = STACK_SIZE;
}

static inline
void voxel_add_atom( voxel * v, atom * a ) {
  if (v->num_atoms==v->max_atoms) {
    v->max_atoms *= 1.3;
    if (v->atoms==v->stack_atoms) {
      v->atoms=(atom**)malloc(v->max_atoms*sizeof(atom *));
      memcpy(v->atoms, v->stack_atoms, STACK_SIZE*sizeof(atom *));
    } else {
      v->atoms=(atom**)realloc(v->atoms, v->max_atoms*sizeof(atom *));
    }
  }
  v->atoms[v->num_atoms++] = a;
}

static inline void voxel_release( voxel * v ) {
  if (v->atoms != v->stack_atoms) free(v->atoms);
}

/* find the bounding box of selected atoms. Return 0 if no atoms 
 * selected. */
static
int find_bbox(const Selection& S, SystemPtr sys, double *min, double *max) {
  double xmin, ymin, zmin;
  double xmax, ymax, zmax;

  // find first selected atom
  int i,n=S.size();
  for (i=0; i<n; i++) if (S[i]) break;
  if (i==n) return 0;

  // initialize min/max to first atom position
  xmin=xmax=sys->atom(i).x;
  ymin=ymax=sys->atom(i).y;
  zmin=zmax=sys->atom(i).z;

  // extend minmax
  for (;i<n; i++) {
    if (!S[i]) continue;
    atom& atm = sys->atom(i);
    double x=atm.x;
    double y=atm.y;
    double z=atm.z;

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
void find_voxel_full_shell_neighbors(voxel *mesh, int nx, int ny, int nz) {
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

    public:
    WithinPredicate( SystemPtr e, double r, bool excl, PredicatePtr s )
      : sys(e), rad(r), sub(s), exclude(excl) {}

    void eval( Selection& s );
    void dump( std::ostream& str ) const {
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

void WithinPredicate::eval( Selection& S ) {

  Selection subsel = full_selection(sys);
  double min[3], max[3];
  double xmin, ymin, zmin;
  double xsize, ysize, zsize;
  const double r2 = rad*rad;
  double ir;
  voxel * mesh;

  int nx, ny, nz, nvoxel;
  int i;

  if (rad <= 0) {
    S.clear();
    return;
  }
  ir = 1.0/rad;

  /* evaluate the subselection */
  sub->eval(subsel);

  /* find bounding box of subselection */
  if (!find_bbox(subsel, sys, min, max)) {
    S.clear();
    return;
  }

  /* If this is an "exwithin" selection, AND S with the converse of the
   * subselection */
  if (exclude) {
    S.subtract(subsel);
  }

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

  mesh = (voxel *)malloc(nvoxel*sizeof(*mesh));
  for (i=0; i<nvoxel; i++) voxel_init(mesh+i);

  find_voxel_full_shell_neighbors(mesh, nx, ny, nz);

  /* map subsel atoms to voxels */
  for (Id i=0; i<subsel.size(); i++) {
    if (!subsel[i]) continue;

    atom& atm = sys->atom(i);
    double x=atm.x;
    double y=atm.y;
    double z=atm.z;
    int xi = (x-xmin)*ir;
    int yi = (y-ymin)*ir;
    int zi = (z-zmin)*ir;
    int index = xi + nx*(yi + ny*zi);
    voxel_add_atom( mesh+index, &atm );
  }

  /* loop over atoms in left selection */
  for (Id i=0; i<S.size(); i++) {
    double x,y,z;
    int xi,yi,zi;
    int j, index;
    int n_nbrs;
    const int * nbrs;
    int on;
    voxel * v;
    if (!S[i]) continue;
    if (subsel[i]) continue;
    atom& atm = sys->atom(i);
    x=atm.x;
    y=atm.y;
    z=atm.z;
    xi = (x-xmin)*ir;
    yi = (y-ymin)*ir;
    zi = (z-zmin)*ir;
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
      voxel * nbr = mesh + nbrs[j];
      atom **atoms = nbr->atoms;
      unsigned k, natoms = nbr->num_atoms;
      for (k=0; k<natoms; k++) {
        double lx=atoms[k]->x;
        double ly=atoms[k]->y;
        double lz=atoms[k]->z;
        double dx = lx-x;
        double dy = ly-y;
        double dz = lz-z;
        if (dx*dx + dy*dy + dz*dz <=r2) {
          on=1;
          break;
        }
      }
      if (on) break;
    }
    if (!on) S[i]=0;
  }
  for (i=0; i<nvoxel; i++) voxel_release(mesh+i);
  free(mesh);
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
    return PredicatePtr(new WithinPredicate(sys,r,false,s));
  }

  PredicatePtr exwithin_predicate( SystemPtr sys, double r, PredicatePtr s ) {
    return PredicatePtr(new WithinPredicate(sys,r,true,s));
  }

  PredicatePtr withinbonds_predicate( SystemPtr sys, int n, PredicatePtr s ) {
    return PredicatePtr(new WithinBondsPredicate(sys,n,s));
  }

}}}

