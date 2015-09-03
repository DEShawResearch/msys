#ifndef desres_pfx_glue_hxx
#define desres_pfx_glue_hxx

#include <algorithm>
#include <vector>
#include <cstdlib>

namespace desres { namespace msys { namespace pfx {

    /* Replace x with number of increments of size 1 or -1
     * to shift each coordinate in order to minimize the square distance
     * between them */
    template <typename scalar>
    void find_optimal_shifts(int n, scalar* x) {
        int i, k;
        scalar* s, *xx, *v, left, tot=0, tot2=0;
        static const scalar MINUS_HALF = -0.5;
        static const scalar PLUS_HALF =  0.5;

        if (n<2) return;
        s = (scalar *)calloc(3*n, sizeof(*s));
        xx = s+n;
        v = xx+n;

        /* bring all points to [-1/2, 1/2). */
        for (i=0; i<n; i++) {
            scalar y = x[i];
            while (y < MINUS_HALF) {
                y += 1;
                s[i] += 1;
            }
            while (y >= PLUS_HALF) {
                y -= 1;
                s[i] -= 1;
            }
            tot += y;
            tot2 += y*y;
            x[i] = xx[i] = y;
        }
        std::sort(xx, xx+n);

        /* We have v(1) = 1/2 sum(xi-xj)^2 = N * sum(xi^2) - sum(xi)^2, 
         * and v(k) = 1/2 sum(xi - xj + (i<k) - (j<k))^2.  The leftmost 
         * point should be the one corresponding to the smallest element 
         * in v.  We can compute v in linear time using a recurrence 
         * relation.  Total complexity O(n log n) for the sort.
         * Credit to R. Lippert for the algorithm.
         */
        v[0] = n*tot2 - tot*tot;
        for (i=0; i<n-1; i++) {
            v[i+1] = v[i] + n-1 - 2*(tot + i - n*xx[i]);
        }
        k = std::min_element(v, v+n)-v;
        left = xx[k];
        for (i=0; i<n; i++) {
            if (x[i] < left) {
                s[i] += 1;
            }
            x[i] = s[i];
        }
        free(s);
    }

    template <typename scalar>
    void apply_glue(unsigned        n,
                    const unsigned *atoms,
                    const unsigned *sizes,
                    const unsigned *ccnts,
                    const scalar   *cell, 
                    const scalar   *proj, 
                    scalar         *pos) {
        if (n<2) return;    // nothing to do
        std::vector<scalar> s(3*n);
        const unsigned *ptr = &atoms[0];
        for (unsigned i=0; i<n; i++) {
            scalar c[3]={0,0,0};
            // compute center using ccnt atoms
            const unsigned sz=ccnts[i];
            scalar ilen = 1/(scalar)sz;
            for (unsigned j=0; j<sz; j++) {
                const scalar* p = pos+3*ptr[j];
                c[0] += p[0];
                c[1] += p[1];
                c[2] += p[2];
            }
            for (unsigned j=0; j<3; j++) c[j] *= ilen;
            // c now holds center of fragment i
            for (unsigned j=0; j<3; j++) {
                const scalar* p = proj+3*j;
                s[i+j*n] = c[0]*p[0] + c[1]*p[1] + c[2]*p[2];
            }
            ptr += sizes[i];
        }

        // s[..n] convert fractional projection to optimal shifts
        for (unsigned i=0; i<3; i++) find_optimal_shifts(n, &s[i*n]);

        // apply shifts to fragments
        ptr = &atoms[0];
        for (unsigned i=0; i<n; i++) {
            scalar sa = s[i];
            scalar sb = s[i+n];
            scalar sc = s[i+2*n];
            scalar dx = cell[0]*sa + cell[3]*sb + cell[6]*sc;
            scalar dy = cell[1]*sa + cell[4]*sb + cell[7]*sc;
            scalar dz = cell[2]*sa + cell[5]*sb + cell[8]*sc;
            const unsigned sz = sizes[i];
            for (unsigned j=0; j<sz; j++) {
                scalar* p = pos+3*ptr[j];
                p[0] += dx;
                p[1] += dy;
                p[2] += dz;
            }
            ptr += sz;
        }
    }

    class Glue {
        // atoms grouped by connected component
        std::vector<unsigned> _atoms;

        // size of each component
        std::vector<unsigned> _sizes;

        // number in each component used for centering
        std::vector<unsigned> _ccnts;

    public:
        typedef unsigned Id;

        Glue(std::vector<Id>& atoms,
             std::vector<Id>& sizes,
             unsigned n,
             const Id* glue) {

            if (n==0 || glue==nullptr || atoms.empty()) return;
            std::vector<Id> clique;
            const Id* ptr;

            // make a lookup table for selected glue atoms
            std::vector<bool> hash(1+*std::max_element(atoms.begin(), atoms.end()));
            for (unsigned i=0; i<n; i++) hash.at(glue[i])=true;

            // find affected components
            ptr = &atoms[0];
            for (unsigned i=0, n=sizes.size(); i<n; i++) {
                unsigned ccnt = 0;
                const unsigned sz = sizes[i];
                for (unsigned j=0; j<sz; j++) {
                    Id atm = ptr[j];
                    if (hash.at(atm)) {
                        _atoms.push_back(atm);
                        ++ccnt;
                    }
                }
                // Found some overlap.  Add the rest of the atoms; i.e. 
                // those that are part of the fragment but not used to 
                // compute the center.
                if (ccnt) {
                    for (unsigned j=0; j<sz; j++) {
                        Id atm = ptr[j];
                        if (!hash[atm]) {
                            _atoms.push_back(atm);
                        }
                    }
                    _sizes.push_back(sz);
                    _ccnts.push_back(ccnt);
                    clique.push_back(i);
                }
                ptr += sz;
            }
            
            // new graph is finished.  Now modify input atoms and sizes.
            if (clique.size()>=2) {
                std::vector<Id> newatoms(_atoms);
                std::vector<Id> newsizes(1, _atoms.size());
                ptr = &atoms[0];
                for (unsigned i=0, n=sizes.size(); i<n; i++) {
                    const unsigned sz = sizes[i];
                    if (!std::binary_search(clique.begin(), clique.end(), i)) {
                        newatoms.insert(newatoms.end(), ptr, ptr+sz);
                        newsizes.push_back(sz);
                    }
                    ptr += sz;
                }
                atoms.swap(newatoms);
                sizes.swap(newsizes);
            }
        }

        // apply unit cell shifts to bring glued fragments together
        template <typename scalar>
        void join(const scalar* cell, const scalar* proj, 
                        scalar* pos) const {
            apply_glue(_sizes.size(),
                       &_atoms[0], &_sizes[0], &_ccnts[0],
                       cell, proj, pos);
        }
    };

}}}

#endif
