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

            template <bool periodic>
            static 
            void find_neighbors(voxel* mesh, int nx, int ny, int nz);
        };

        void rotate(const Float* mat, const Float* p, Float* q) {
            Float x=p[0];
            Float y=p[1];
            Float z=p[2];
            q[0]=x*mat[0] + y*mat[1] + z*mat[2];
            q[1]=x*mat[3] + y*mat[4] + z*mat[5];
            q[2]=x*mat[6] + y*mat[7] + z*mat[8];
        }
    }

    template <bool pbc>
    void details::voxel::find_neighbors(
            voxel* mesh, int nx, int ny, int nz) {
        int r,s,t;
        for (int zi=0; zi<nz; zi++) {
            for (int yi=0; yi<ny; yi++) {
                for (int xi=0; xi<nx; xi++) {
                    int self=xi + nx*(yi+ny*zi);
                    int* nbrs = mesh[self].nbrs;
                    int n=0;
                    for (int ti=zi-1; ti<=zi+1; ti++) {
                        if (pbc) t = ti<0 ? nz-1 : ti==nz ? 0 : ti;
                        else { if (ti<0 || ti>=nz) continue; t=ti; }
                        for (int si=yi-1; si<=yi+1; si++) {
                            if (pbc) s = si<0 ? ny-1 : si==ny ? 0 : si;
                            else { if (si<0 || si>=ny) continue; s=si; }
                            for (int ri=xi-1; ri<=xi+1; ri++) {
                                if (pbc) r = ri<0 ? nx-1 : ri==nx ? 0 : ri;
                                else { if (ri<0 || ri>=nx) continue; r=ri; }
                                int index = r+nx*(s+ny*t);
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
                       const Float* cell,
                       Iter beginA, Iter endA,
                       Iter beginB, Iter endB,
                       Output const& output) {

        using namespace details;

        if (rad<=0 || beginA==endA || beginB==endB || !pos) return;
        Float min[3], max[3], rot[9];
        std::vector<Float> rotated_positions;

        if (cell) {
            /* If non-orthorhombic cell, make a copy of the positions and
             * rotate them to align them with the coordinate axes.  */
            bool orthorhombic = true;
            for (int i=0; i<3; i++) for (int j=0; j<3; j++) {
                if (cell[3*i+j]!=0 && i!=j) {
                    orthorhombic = false;
                    break;
                }
            }
            Float r[3];
            /* generate rotation matrix */
            for (int i=0; i<3; i++) {
                Float x=cell[3*i+0];
                Float y=cell[3*i+1];
                Float z=cell[3*i+2];
                r[i] = sqrt(x*x + y*y + z*z);
                if (r[i]>0) {
                    Float ir=1/r[i];
                    x *= ir;
                    y *= ir;
                    z *= ir;
                }
                rot[3*i+0]=x;
                rot[3*i+1]=y;
                rot[3*i+2]=z;
            }

            /* place atoms inside rotated cell */
            IdList ids;
            ids.insert(ids.end(), beginA, endA);
            ids.insert(ids.end(), beginB, endB);
            for (IdList::const_iterator iter=ids.begin(), e=ids.end(); iter!=e; ++iter) {
                Id atm = *iter;
                if (rotated_positions.size()<=3*atm) {
                    rotated_positions.resize(3*atm+3);
                }
                if (orthorhombic) {
                    memcpy(&rotated_positions[3*atm], pos+3*atm,
                            3*sizeof(Float));
                } else {
                    rotate(rot, pos+3*atm, &rotated_positions[3*atm]);
                }

                for (int i=0; i<3; i++) {
                    Float c = rotated_positions[3*atm+i];
                    Float m = r[i];
                    while (c<0) c+=m;
                    while (c>=r[i]) c-=m;
                    rotated_positions[3*atm+i]=c;
                    /* bounding box from cell dimensions */
                    min[i] = 0;
                    max[i] = m;
                }
            }
            pos = &rotated_positions[0];

        } else {
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
        if (cell) {
            voxel::find_neighbors<true>(mesh, nx, ny, nz);
        } else {
            voxel::find_neighbors<false>(mesh, nx, ny, nz);
        }

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
            if (!cell && (xi<0 || xi>=nx || yi<0 || yi>=ny || zi<0 || zi>=nz)) {
                continue;
            }
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
                    if (cell) {
                        /* periodic distance from q to p */
                        for (int i=0; i<3; i++) {
                            Float d = p[i]-q[i];
                            Float m = max[i];
                            Float m2 = m/2;
                            while (d>m2) d-=m;
                            while (d<-m2) d+=m;
                            d2 += d*d;
                        }

                    } else {
                        Float dx=x-q[0];
                        Float dy=y-q[1];
                        Float dz=z-q[2];
                        d2 = dx*dx + dy*dy + dz*dz;
                    }
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
