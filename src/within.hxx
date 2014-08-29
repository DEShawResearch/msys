#ifndef desres_msys_within_hxx
#define desres_msys_within_hxx

#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <vector>
#include <limits>
#include <algorithm>

#ifdef __SSE2__
#include <emmintrin.h>
#endif

#ifdef __SSE4_1__
#include <smmintrin.h>
#endif

static inline __m128i quad_int_multiply(const __m128i &a, const __m128i &b) {
#ifdef __SSE4_1__
    return _mm_mullo_epi32(a, b);
#else
    __m128i tmp1 = _mm_mul_epu32(a,b); /* mul 2,0*/
    __m128i tmp2 = _mm_mul_epu32( _mm_srli_si128(a,4), _mm_srli_si128(b,4)); /* mul 3,1 */
    return _mm_unpacklo_epi32(_mm_shuffle_epi32(tmp1, _MM_SHUFFLE (0,0,2,0)), _mm_shuffle_epi32(tmp2, _MM_SHUFFLE (0,0,2,0))); /* shuffle results to [63..0] and pack */
#endif
}

#include "types.hxx"

/*
 * Voxel grid:
 *
 * Number of voxels in each dimension sufficient to bound the input data
 * is given by normalized xmax divided by adjusted voxel size, rounded up.
 *
 * I create a dummy voxel at the left and right of each dimension.
 *
 * I also want to create extra space for voxels so that I don't have to 
 * worry about handling the edge voxels differently.  That works out to
 * be an extra slab of ny*nz voxels for the negative x coordinate, an 
 * extra row of nz voxels for the negative y coordinate, and an extra single
 * voxel for the negative z.  Same thing at the other end.  I don't have to
 * worry about the internal edge voxels.  If I'm at x, y, or z=0 and I access
 * the voxel to the "left", I'll just hit a dummy voxel.
 *
 */

namespace desres { namespace msys {

    template <typename scalar>
    class SpatialHash {

        /* Support distance queries of at most _radius */
        scalar rad;
        
        /* inverse radius */
        scalar ir;

        /* bounding box of target coordinates */
        scalar xmin, ymin, zmin;
        scalar xmax, ymax, zmax;

        /* location of grid origin */
        scalar ox, oy, oz;

        /* voxel grid dimension */
        int nx, ny, nz;

        /* neighbor full shell offsets */
        int full_shell[10];
        int strip_lens[10];

        /* target coordinates */
        const int ntarget;
        scalar *_x, *_y, *_z;
        scalar *_tmpx, *_tmpy, *_tmpz;

        /* cumulative number of particles hashed to the i'th voxel */
        std::vector<uint32_t> _countspace;
        uint32_t *_counts;

        void compute_full_shell() {
            int* p = full_shell;
            int* n = strip_lens;
            /* self and z+1 */
            *p++ = 0;   /* self */
            *n++ = 2;
            /* z-1 */
            *p++ = -1;
            *n++ = 1;
            /* all others as groups of three strips */
            for (int i=-1; i<=1; i++) {
                for (int j=-1; j<=1; j++) {
                    if (i!=0 || j!=0) {
                        *p++ = -1 + nz*(j + ny*i);
                        *n++ = 3;
                    }
                }
            }
        }


    public:
        ~SpatialHash() {
            free(_x);
            free(_y);
            free(_z);
            free(_tmpx);
            free(_tmpy);
            free(_tmpz);
        }

        SpatialHash( const scalar *pos, IdList const& ids)
        : rad(), ir(), ntarget(ids.size()), _x(), _y(), _z(), _counts() {
            nx = ny = nz = 0;
            ox = oy = oz = 0;
            if (ntarget<1) return;

            /* copy to transposed arrays */
            posix_memalign((void **)&_x, 16, ntarget*sizeof(*_x));
            posix_memalign((void **)&_y, 16, ntarget*sizeof(*_y));
            posix_memalign((void **)&_z, 16, ntarget*sizeof(*_z));
            posix_memalign((void **)&_tmpx, 16, ntarget*sizeof(*_tmpx));
            posix_memalign((void **)&_tmpy, 16, ntarget*sizeof(*_tmpy));
            posix_memalign((void **)&_tmpz, 16, ntarget*sizeof(*_tmpz));
            for (int i=0; i<ntarget; i++) {
                const scalar *xyz = pos+3*ids[i];
                _x[i] = xyz[0];
                _y[i] = xyz[1];
                _z[i] = xyz[2];
            }

            /* compute bounds for positions */
            find_bbox(ntarget, _x, &xmin, &xmax);
            find_bbox(ntarget, _y, &ymin, &ymax);
            find_bbox(ntarget, _z, &zmin, &zmax);
        }

        static
        void find_bbox(int n, const scalar* x, scalar *_min, scalar *_max) {
            scalar min = x[0], max=x[0];
            for (int i=0; i<n; i++) {
                scalar p = x[i];
                min = min < p ? min : p;
                max = max > p ? max : p;
            }
            *_min = min;
            *_max = max;
        }

        SpatialHash& voxelize(scalar r) {
            rad = r;
            ir = scalar(1)/rad;
            ox = xmin - rad;
            oy = ymin - rad;
            oz = zmin - rad;

            /* construct voxel grid.  */
            nx = (xmax-xmin)*ir + 3;
            ny = (ymax-ymin)*ir + 3;
            nz = (zmax-zmin)*ir + 3;
            int nvoxels = nx*ny*nz;

            /* allocate space for extra voxels so we don't have to worry
             * about edge cases. */
            int voxpad = 1 + nz*(1 + ny);
            _countspace.resize(nvoxels+2*voxpad+1);
            std::fill(_countspace.begin(), _countspace.end(), 0);
            _counts = &_countspace[voxpad];
            ++_counts;
            compute_full_shell();

            /* map points to voxels and compute voxel histogram */
            /* FIXME: SIMD */
            std::vector<uint32_t> voxids(ntarget);
            for (int i=0; i<ntarget; i++) {
                scalar x = _x[i];
                scalar y = _y[i];
                scalar z = _z[i];
                int xi = (x-ox) * ir;
                int yi = (y-oy) * ir;
                int zi = (z-oz) * ir;
                int voxid = zi + nz*(yi + ny*xi);
                voxids[i] = voxid;
                ++_counts[voxid];
            }

            /* compute starting index for each voxel */
            uint32_t cnt = 0;
            for (int i=0; i<=nvoxels; i++) {
                uint32_t tmp = _counts[i];
                _counts[i] = cnt;
                cnt += tmp;
            }

            /* copy positions sorted by voxel */
            for (int i=0; i<ntarget; i++) {
                int voxid = voxids[i];
                unsigned& j = _counts[voxid];
                _tmpx[j] = _x[i];
                _tmpy[j] = _y[i];
                _tmpz[j] = _z[i];
                ++j;
            }
            std::swap(_x,_tmpx);
            std::swap(_y,_tmpy);
            std::swap(_z,_tmpz);

            /* shift counts up by to undo the counting sort we just did */
            --_counts;
            return *this;
        }

        /* Is the given point within r of any point in the region? */
        bool test(scalar r, scalar x, scalar y, scalar z) const {
            int xi = (x-ox) * ir;
            int yi = (y-oy) * ir;
            int zi = (z-oz) * ir;
            if (xi<0 || xi>=nx ||
                yi<0 || yi>=ny ||
                zi<0 || zi>=nz) return false;
            int voxid = zi + nz*(yi + ny*xi);
            return test2(r*r, voxid, x,y,z);
        }

        scalar mindist2(scalar x, scalar y, scalar z) const {
            int xi = (x-ox) * ir;
            int yi = (y-oy) * ir;
            int zi = (z-oz) * ir;
            scalar r2 = std::numeric_limits<scalar>::max();
            if (xi<0 || xi>=nx ||
                yi<0 || yi>=ny ||
                zi<0 || zi>=nz) return r2;
            int voxid = zi + nz*(yi + ny*xi);
            /* TODO: SIMD */
            for (int i=0; i<10; i++) {
                int vox = voxid + full_shell[i];
                uint32_t b = _counts[vox], e = _counts[vox+strip_lens[i]];
                const scalar *xi = _x+b;
                const scalar *yi = _y+b;
                const scalar *zi = _z+b;

                for (; b<e; ++b, ++xi, ++yi, ++zi) {
                    scalar dx = x - *xi;
                    scalar dy = y - *yi;
                    scalar dz = z - *zi;
                    scalar d2 = dx*dx + dy*dy + dz*dz;
                    r2 = std::min(r2, d2);
                }
            }
            return r2;
        }

        IdList findNearest(unsigned k, const scalar* pos, IdList const& ids,
                          const double* cell = NULL) {

            /* enough atoms available */
            if (ids.size() <= k) return ids;

            IdList smin, smax;
            scalar rmin = 0.0;
            scalar rmax = 2.5;

            /* increase rmax until we get at least K atoms */
            for (;;) {
                smax = findWithin(rmax, pos, ids, cell);
                if (smax.size() >= k) break;
                smin.swap(smax);
                rmin = rmax;
                rmax *= 1.5;
            }

            /* Refine with a couple rounds of bisection search */
            for (int nb=0; nb<6; nb++) {
                scalar rm = 0.5*(rmin+rmax);
                IdList sm = find_within(rm, pos, smax, cell);
                if (sm.size() >= k) {
                    smax.swap(sm);
                    rmax = rm;
                } else {
                    smin.swap(sm);
                    rmin = rm;
                }
            }

            /* get the closest points in smax but not in smin */
            std::vector<std::pair<scalar, Id> > pts;
            pts.reserve(smax.size()-smin.size());
            for (unsigned i=0, n=smax.size(); i<n; i++) {
                Id id = smax[i];
                if (std::binary_search(smin.begin(), smin.end(), id)) continue;
                const scalar* p = pos+3*id;
                scalar r2 = mindist2(p[0], p[1], p[2]);
                pts.push_back(std::make_pair(r2, id));
            }
            unsigned n = k-smin.size();
            std::partial_sort(pts.begin(), pts.begin()+n, pts.end());

            /* append these to smin and return */
            for (unsigned i=0; i<n; i++) smin.push_back(pts[i].second);
            std::sort(smin.begin(), smin.end());
            return smin;
        }
        
        /* Find which of n points are within r of the hashed target set.
         * If ids is NULL, read n contiguous xyz positions from pos and 
         * do not store results.  Otherwise, read the i'th position from 
         * pos + 3*ids[i], and on return store into the start of ids the
         * points which satisfy the distance check.  Return the number found.
         */
        IdList findWithin(scalar r, const scalar* pos, IdList const& ids,
                          const double* cell = NULL) {
            voxelize(r);
            return find_within(r,pos,ids,cell);
        }

        IdList find_within(scalar r, const scalar* pos, IdList const& ids,
                          const double* cell ) const {
            IdList result;
            result.reserve(ids.size());
            const scalar cx = cell ? cell[0] : 0;
            const scalar cy = cell ? cell[4] : 0;
            const scalar cz = cell ? cell[8] : 0;
            int j=0, n=ids.size();
#ifdef __SSE2__
            const scalar r2 = r*r;
            int b0, b1, b2, b3;
            __m128 xm = _mm_set1_ps(ox);
            __m128 ym = _mm_set1_ps(oy);
            __m128 zm = _mm_set1_ps(oz);
            __m128 ir = _mm_set1_ps(this->ir);
           
            __m128i nx  = _mm_set1_epi32(this->nx);
            __m128i ny  = _mm_set1_epi32(this->ny);
            __m128i nz  = _mm_set1_epi32(this->nz);
            __m128i z1  = _mm_set1_epi32(-1);
            __m128i nyz = _mm_set1_epi32(this->nz * this->ny);

            for (; j<n-4; j+=4) {
                __m128 p0, p1, p2;
                __m128 q0, q1, q2;
        
                __m128i i0, i1, i2;
                __m128i j0, j1, j2;
                __m128i c0, c1, c2;

                const scalar *a = pos+3*ids[j  ];
                const scalar *b = pos+3*ids[j+1];
                const scalar *c = pos+3*ids[j+2];
                const scalar *d = pos+3*ids[j+3];
                p0 = _mm_setr_ps(a[0], b[0], c[0], d[0]);
                p1 = _mm_setr_ps(a[1], b[1], c[1], d[1]);
                p2 = _mm_setr_ps(a[2], b[2], c[2], d[2]);

                q0 = _mm_mul_ps(_mm_sub_ps(p0, xm), ir);    /* (x-p.xm)*ir */
                q1 = _mm_mul_ps(_mm_sub_ps(p1, ym), ir);
                q2 = _mm_mul_ps(_mm_sub_ps(p2, zm), ir);
        
                i0 = _mm_cvttps_epi32(q0); /* truncate */
                i1 = _mm_cvttps_epi32(q1);
                i2 = _mm_cvttps_epi32(q2);
        
                j2 = quad_int_multiply(i0, nyz);
                j1 = quad_int_multiply(i1, nz);
        
                j0 = _mm_add_epi32(i2, j1);
                j0 = _mm_add_epi32(j0, j2);
        
                c0 = _mm_and_si128(
                        _mm_cmplt_epi32(z1, i0),    /* -1 < i0 */
                        _mm_cmplt_epi32(i0, nx));   /* i0 < nx */
        
                c1 = _mm_and_si128(
                        _mm_cmplt_epi32(z1, i1),   /* i1 >= 0 */
                        _mm_cmplt_epi32(i1, ny));
        
                c2 = _mm_and_si128(
                        _mm_cmplt_epi32(z1, i2),   /* i2 >= 0 */
                        _mm_cmplt_epi32(i2, nz));
        
                c1 = _mm_and_si128(c0, c1);
                c2 = _mm_and_si128(c2, c1);
        
                j0 = _mm_or_si128(
                        _mm_and_si128   (c2, j0),
                        _mm_andnot_si128(c2, z1));

                if (_mm_movemask_epi8(_mm_cmpgt_epi32(j0, z1))==0) continue;

#ifdef __SSE4_1__
                int v0 = _mm_extract_epi32(j0,0);
                int v1 = _mm_extract_epi32(j0,1);
                int v2 = _mm_extract_epi32(j0,2);
                int v3 = _mm_extract_epi32(j0,3);
#else
                int v0 = _mm_cvtsi128_si32(j0);
                int v1 = _mm_cvtsi128_si32(_mm_shuffle_epi32(j0,0x55));
                int v2 = _mm_cvtsi128_si32(_mm_shuffle_epi32(j0,0xAA));
                int v3 = _mm_cvtsi128_si32(_mm_shuffle_epi32(j0,0xFF));
#endif

#define GETX(i) _mm_cvtss_f32(_mm_shuffle_ps(p0,p0, _MM_SHUFFLE(0, 0, 0, i)))
#define GETY(i) _mm_cvtss_f32(_mm_shuffle_ps(p1,p1, _MM_SHUFFLE(0, 0, 0, i)))
#define GETZ(i) _mm_cvtss_f32(_mm_shuffle_ps(p2,p2, _MM_SHUFFLE(0, 0, 0, i)))
                b0 = (v0>=0) && test2(r2, v0, GETX(0), GETY(0), GETZ(0));
                b1 = (v1>=0) && test2(r2, v1, GETX(1), GETY(1), GETZ(1));
                b2 = (v2>=0) && test2(r2, v2, GETX(2), GETY(2), GETZ(2));
                b3 = (v3>=0) && test2(r2, v3, GETX(3), GETY(3), GETZ(3));
                if (cell) {
                    if (!b0) b0 = minimage(r,cx,cy,cz, GETX(0),GETY(0),GETZ(0));
                    if (!b1) b1 = minimage(r,cx,cy,cz, GETX(1),GETY(1),GETZ(1));
                    if (!b2) b2 = minimage(r,cx,cy,cz, GETX(2),GETY(2),GETZ(2));
                    if (!b3) b3 = minimage(r,cx,cy,cz, GETX(3),GETY(3),GETZ(3));
                }
#undef GETX
#undef GETY
#undef GETZ
                if (b0) result.push_back(ids[j  ]);
                if (b1) result.push_back(ids[j+1]);
                if (b2) result.push_back(ids[j+2]);
                if (b3) result.push_back(ids[j+3]);
            }
#endif
            for (; j<n; j++) {
                unsigned id = ids[j];
                const scalar *xyz = pos + 3*id;
                if (             test(r, xyz[0], xyz[1], xyz[2]) ||
                    minimage(r,cx,cy,cz, xyz[0], xyz[1], xyz[2])) 
                    result.push_back(id);
            }
            return result;
        }

        bool minimage(scalar r, scalar ga, scalar gb, scalar gc,
                                scalar px, scalar py, scalar pz) const {
            scalar xlo = xmin - r;
            scalar ylo = ymin - r;
            scalar zlo = zmin - r;
            scalar xhi = xmax + r;
            scalar yhi = ymax + r;
            scalar zhi = zmax + r;
            for (int i=-1; i<=1; i++) {
                scalar x = px + ga*i;
                if (x<xlo || x>xhi) continue;
                for (int j=-1; j<=1; j++) {
                    scalar y = py + gb*j;
                    if (y<ylo || y>yhi) continue;
                    for (int k=-1; k<=1; k++) {
                        scalar z = pz + gc*k;
                        if (z<zlo || z>zhi) continue;
                        if (i==0 && j==0 && k==0) continue;
                        if (test(r,x,y,z)) return true;
                    }
                }
            }
            return false;
        }


        bool test2(scalar r2, int voxid, scalar x, scalar y, scalar z) const {
#ifdef __SSE2__
            __m128 xj = _mm_set1_ps(x);
            __m128 yj = _mm_set1_ps(y);
            __m128 zj = _mm_set1_ps(z);
            __m128 R2 = _mm_set1_ps(r2);
#endif
            for (int i=0; i<10; i++) {
                int vox = voxid + full_shell[i];
                uint32_t b = _counts[vox], e = _counts[vox+strip_lens[i]];
                const scalar * xi = _x+b;
                const scalar * yi = _y+b;
                const scalar * zi = _z+b;

#ifdef __SSE2__
                /* advance to aligned offset */
                for (; b<e && (b&3); ++b, ++xi, ++yi, ++zi) {
                    scalar dx = x - *xi;
                    scalar dy = y - *yi;
                    scalar dz = z - *zi;
                    scalar d2 = dx*dx + dy*dy + dz*dz;
                    if (d2<=r2) return true;
                }

                /* simd for chunks of four */
                for (; b+4<e; b+=4, xi+=4, yi+=4, zi+=4) {
                    __m128 p0, p1, p2;
                    __m128 q0, q1, q2;

                    p0 = _mm_load_ps(xi);
                    p1 = _mm_load_ps(yi);
                    p2 = _mm_load_ps(zi);
        
                    q0 = _mm_sub_ps(p0, xj);    /* dx */
                    q1 = _mm_sub_ps(p1, yj);    /* dy */
                    q2 = _mm_sub_ps(p2, zj);    /* dz */
                    
                    p0 = _mm_mul_ps(q0, q0);    /* dx**2 */
                    p1 = _mm_mul_ps(q1, q1);    /* dy**2 */
                    p2 = _mm_mul_ps(q2, q2);    /* dz**2 */

                    q0 = _mm_add_ps(_mm_add_ps(p0,p1),p2);
                    if (_mm_movemask_ps(_mm_cmple_ps(q0, R2))) return true;
                }
#endif

                /* stragglers */
                for (; b<e; ++b, ++xi, ++yi, ++zi) {
                    scalar dx = x - *xi;
                    scalar dy = y - *yi;
                    scalar dz = z - *zi;
                    scalar d2 = dx*dx + dy*dy + dz*dz;
                    if (d2<=r2) return true;
                }
            }
            return false;
        }

        scalar radius() const { return rad; }
    };
}}

#endif

