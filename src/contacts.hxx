#ifndef desres_msys_contacts_hxx
#define desres_msys_contacts_hxx

#include "types.hxx"
#include <string.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>

namespace desres { namespace msys {

    namespace details {

        struct voxel {
            static const int INIT_STACK_SIZE = 4;
            Id* pts;
            int num;
            int max;
            int n_nbrs;
            int nbrs[27];
            Id stack[INIT_STACK_SIZE];

            voxel() : pts(stack), num(0), max(INIT_STACK_SIZE), n_nbrs(0) {}
            ~voxel() { if (pts!=stack) free(pts); }
            void add(Id pt) {
                if (num==max) {
                    max *= 1.34;
                    if (pts==stack) {
                        pts = (Id *)malloc(max*sizeof(*pts));
                        memcpy(pts, stack, num*sizeof(*pts));
                    } else {
                        pts = (Id *)realloc(pts, max*sizeof(*pts));
                    }
                }
                pts[num++] = pt;
            }

            static 
            void find_neighbors(voxel* mesh, int nx, int ny, int nz);
        };

    }

    void details::voxel::find_neighbors(
            voxel* mesh, int nx, int ny, int nz) {
        for (int zi=0; zi<nz; zi++) {
            for (int yi=0; yi<ny; yi++) {
                for (int xi=0; xi<nx; xi++) {
                    int self=xi + nx*(yi+ny*zi);
                    int* nbrs = mesh[self].nbrs;
                    int n=0;
                    for (int ti=zi-1; ti<=zi+1; ti++) {
                        if (ti<0 || ti>=nz) continue;
                        for (int si=yi-1; si<=yi+1; si++) {
                            if (si<0 || si>=ny) continue;
                            for (int ri=xi-1; ri<=xi+1; ri++) {
                                if (ri<0 || ri>=nx) continue;
                                int index = ri+nx*(si+ny*ti);
                                if (mesh[index].num) nbrs[n++] = index;
                            }
                        }
                    }
                    mesh[self].n_nbrs = n;
                }
            }
        }
    }

    template <typename Float, typename Iter, typename Output>
    void find_contacts(Float rad, 
                       const Float* pos,
                       Iter beginA, Iter endA,
                       Iter beginB, Iter endB,
                       Output const& output) {

        using namespace details;

        if (rad<=0 || beginA==endA || beginB==endB || !pos) return;
        Float min[3], max[3];

        /* bounding box from subselection */
        const Float* p = pos + *beginB;
        for (int i=0; i<3; i++) min[i] = max[i] = p[i];
        for (Iter iter=beginB, e=endB; iter!=e; ++iter) {
            p = pos + 3*(*iter);
            for (int i=0; i<3; i++) {
                Float c = p[i];
                min[i] = std::min(min[i], c);
                max[i] = std::max(max[i], c);
            }
        }
        /* extend bounding box by selection radius */
        for (int i=0; i<3; i++) {
            min[i] -= rad;
            max[i] += rad;
        }

        /* construct voxel mesh */
        Float xmin = min[0];
        Float ymin = min[1];
        Float zmin = min[2];
        Float xsize = max[0]-xmin;
        Float ysize = max[1]-ymin;
        Float zsize = max[2]-zmin;
        const Float ir = 1/rad;
        const Float r2 = rad*rad;
        int nx = (int)(xsize*ir)+1;
        int ny = (int)(ysize*ir)+1;
        int nz = (int)(zsize*ir)+1;
        int nvoxel = nx*ny*nz;

        voxel* mesh = new voxel[nvoxel];

        /* map B atoms to voxels */
        for (Iter iter=beginB, e=endB; iter!=e; ++iter) {
            Id atm = *iter;
            const Float* p = pos+3*atm;
            int xi = (p[0]-xmin)*ir;
            int yi = (p[1]-ymin)*ir;
            int zi = (p[2]-zmin)*ir;
            int index = xi + nx*(yi+ny*zi);
            assert(index>=0 && index<=nvoxel);
            mesh[index].add(atm);
        }

        /* set up periodicity in mesh */
        voxel::find_neighbors(mesh, nx, ny, nz);

        /* find contacts with A atoms */
        for (Iter iter=beginA, e=endA; iter!=e; ++iter) {
            Id atm = *iter;
            const Float* p = pos+3*atm;
            Float x = p[0];
            Float y = p[1];
            Float z = p[2];
            int xi = (x-xmin)*ir;
            int yi = (y-ymin)*ir;
            int zi = (z-zmin)*ir;
            if (xi<0 || xi>=nx || yi<0 || yi>=ny || zi<0 || zi>=nz) continue;
            const int index = xi + nx*(yi + ny*zi);
            const voxel& v = mesh[index];
            const int n_nbrs = v.n_nbrs;
            const int* nbrs = v.nbrs;
            for (int j=0; j<n_nbrs; j++) {
                const voxel& nbr = mesh[nbrs[j]];
                const int natoms = nbr.num;
                for (int k=0; k<natoms; k++) {
                    const Id pk = nbr.pts[k];
                    if (atm==pk || output.exclude(atm, pk)) continue;
                    const Float* q = pos+3*pk;
                    Float d2=0;
                    Float dx=x-q[0];
                    Float dy=y-q[1];
                    Float dz=z-q[2];
                    d2 = dx*dx + dy*dy + dz*dz;
                    if (d2 <= r2) {
                        output(atm, pk, d2);
                    }
                }
            }
        }

        delete[] mesh;

    }




}}

#endif
