#include "spatial_hash.hxx"
#include <limits>
#include <math.h>
#include "pfx/rms.hxx"
#include <stdio.h>

#define EXCLUDE_SELF_CONTACTS

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

static
void find_bbox(int n, const float* x, float *_min, float *_max) {
    float min = x[0], max=x[0];
    for (int i=0; i<n; i++) {
        float p = x[i];
        min = min < p ? min : p;
        max = max > p ? max : p;
    }
    *_min = min;
    *_max = max;
}

using namespace desres::msys;

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

void SpatialHash::compute_full_shell() {
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

SpatialHash::~SpatialHash() {
    free(_x);
    free(_y);
    free(_z);
    free(_tmpx);
    free(_tmpy);
    free(_tmpz);
    free(rot);
    free(_ids);
    free(_tmpids);
}

SpatialHash::SpatialHash( const float *pos, int n, const Id* ids, 
                          const double* cell)
: rad(), ir(), 
  xmin(), ymin(), zmin(),
  xmax(), ymax(), zmax(),
  rot(), cx(), cy(), cz(),
  ntarget(n), 
  _x(), _y(), _z(), 
  _tmpx(), _tmpy(), _tmpz(), 
  _counts(), _ids(), _tmpids() {

    nx = ny = nz = 0;
    ox = oy = oz = 0;
    if (ntarget<1) return;

    if (cell) {
        cx = sqrt(cell[0]*cell[0] + cell[1]*cell[1] + cell[2]*cell[2]);
        cy = sqrt(cell[3]*cell[3] + cell[4]*cell[4] + cell[5]*cell[5]);
        cz = sqrt(cell[6]*cell[6] + cell[7]*cell[7] + cell[8]*cell[8]);
        if (cx==0 || cy==0 || cz==0) {
            MSYS_FAIL("cell has zero-length dimensions");
        }
        posix_memalign((void **)&rot, 16, 9*sizeof(*rot));
        double d1=0, d2=0, d3=0; /* row dot-products */
        for (int i=0; i<3; i++) {
            rot[0+i] = cell[0+i]/cx;
            rot[3+i] = cell[3+i]/cy;
            rot[6+i] = cell[6+i]/cz;
            d1 += rot[0+i]*rot[3+i];
            d2 += rot[0+i]*rot[6+i];
            d3 += rot[3+i]*rot[6+i];
        }
        static const float eps = 1e-4;
        if (fabs(d1)>eps || fabs(d2)>eps || fabs(d3)>eps) {
            MSYS_FAIL("cell appears triclinic: dot products " << d1 << " " << d2 << " " << d3);
        }
        if (!(rot[1] || rot[2] || rot[3] || rot[5] || rot[6] || rot[7])) {
            free(rot);
            rot=NULL;
        }
    }

    /* copy to transposed arrays */
    posix_memalign((void **)&_x, 16, ntarget*sizeof(*_x));
    posix_memalign((void **)&_y, 16, ntarget*sizeof(*_y));
    posix_memalign((void **)&_z, 16, ntarget*sizeof(*_z));
    posix_memalign((void **)&_tmpx, 16, ntarget*sizeof(*_tmpx));
    posix_memalign((void **)&_tmpy, 16, ntarget*sizeof(*_tmpy));
    posix_memalign((void **)&_tmpz, 16, ntarget*sizeof(*_tmpz));
    posix_memalign((void **)&_ids, 16, ntarget*sizeof(*_ids));
    posix_memalign((void **)&_tmpids, 16, ntarget*sizeof(*_tmpids));
    for (int i=0; i<ntarget; i++) {
        Id id = ids ? ids[i] : i;
        const float *xyz = pos+3*id;
        if(rot) {
            float pos[3]={xyz[0],xyz[1],xyz[2]};
            pfx::apply_rotation(1,pos,rot);
            _x[i] = pos[0];
            _y[i] = pos[1];
            _z[i] = pos[2];
        } else {
            _x[i] = xyz[0];
            _y[i] = xyz[1];
            _z[i] = xyz[2];
        }
        _ids[i] = id;
    }

    /* compute bounds for positions */
    find_bbox(ntarget, _x, &xmin, &xmax);
    find_bbox(ntarget, _y, &ymin, &ymax);
    find_bbox(ntarget, _z, &zmin, &zmax);
}


SpatialHash& SpatialHash::voxelize(float r) {
    if (r<=0) MSYS_FAIL("radius " << r << " must be positive");
    rad = r;
    ir = float(1)/rad;
    ox = xmin - rad;
    oy = ymin - rad;
    oz = zmin - rad;

    /* construct voxel grid.  */
    nx = (xmax-xmin)*ir + 3;
    ny = (ymax-ymin)*ir + 3;
    nz = (zmax-zmin)*ir + 3;
    static const int maxdim = 500;
    if (nx > maxdim || ny > maxdim || nz > maxdim) {
        float dbig = std::max(std::max(xmax-xmin, ymax-ymin), zmax-zmin);
        ir = float(maxdim) / dbig;
        nx = (xmax-xmin)*ir + 3;
        ny = (ymax-ymin)*ir + 3;
        nz = (zmax-zmin)*ir + 3;
    }
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
        float x = _x[i];
        float y = _y[i];
        float z = _z[i];
        int xi = (x-ox) * ir;
        int yi = (y-oy) * ir;
        int zi = (z-oz) * ir;
        int voxid = zi + nz*(yi + ny*xi);
        voxids[i] = voxid;
        ++_counts[voxid];
    }

    /* compute starting index for each voxel */
    uint32_t cnt = 0;
    maxcount = 0;
    for (int i=0; i<=nvoxels; i++) {
        uint32_t tmp = _counts[i];
        maxcount = std::max(maxcount, tmp);
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
        _tmpids[j] = _ids[i];
        ++j;
    }
    std::swap(_x,_tmpx);
    std::swap(_y,_tmpy);
    std::swap(_z,_tmpz);
    std::swap(_ids, _tmpids);

    /* space for contacts */
    maxcount += 3;  // since we're writing in chunks of four

    /* shift counts up by to undo the counting sort we just did */
    --_counts;
    return *this;
}

float SpatialHash::mindist2(float x, float y, float z) const {
    int xi = (x-ox) * ir;
    int yi = (y-oy) * ir;
    int zi = (z-oz) * ir;
    float r2 = std::numeric_limits<float>::max();
    if (xi<0 || xi>=nx ||
        yi<0 || yi>=ny ||
        zi<0 || zi>=nz) return r2;
    int voxid = zi + nz*(yi + ny*xi);
    /* TODO: SIMD */
    for (int i=0; i<10; i++) {
        int vox = voxid + full_shell[i];
        uint32_t b = _counts[vox], e = _counts[vox+strip_lens[i]];
        const float *xi = _x+b;
        const float *yi = _y+b;
        const float *zi = _z+b;

        for (; b<e; ++b, ++xi, ++yi, ++zi) {
            float dx = x - *xi;
            float dy = y - *yi;
            float dz = z - *zi;
            float d2 = dx*dx + dy*dy + dz*dz;
            r2 = std::min(r2, d2);
        }
    }
    return r2;
}

IdList SpatialHash::findNearest(unsigned k, const float* pos, 
                   unsigned num, const Id* ids) {

    /* If there are no hashed atoms, we should stop now. */
    if (ntarget==0) MSYS_FAIL("No atoms in target selection");

    /* enough atoms available */
    if (num <= k) return IdList(ids, ids+num);

    IdList smin, smax;
    float rmin = 0.0;
    float rmax = 2.5;

    /* increase rmax until we get at least K atoms */
    for (;;) {
        smax = findWithin(rmax, pos, num, ids);
        if (smax.size() >= k) break;
        smin.swap(smax);
        rmin = rmax;
        rmax *= 1.5;
    }

    /* Refine with a couple rounds of bisection search */
    for (int nb=0; nb<6; nb++) {
        float rm = 0.5*(rmin+rmax);
        IdList sm = find_within(rm, pos, smax.size(),&smax[0]);
        if (sm.size() >= k) {
            smax.swap(sm);
            rmax = rm;
        } else {
            smin.swap(sm);
            rmin = rm;
        }
    }

    /* get the closest points in smax but not in smin */
    std::vector<std::pair<float, Id> > pts;
    pts.reserve(smax.size()-smin.size());
    for (unsigned i=0, n=smax.size(); i<n; i++) {
        Id id = smax[i];
        if (std::binary_search(smin.begin(), smin.end(), id)) continue;
        const float* p = pos+3*id;
        float r2 = mindist2(p[0], p[1], p[2]);
        pts.push_back(std::make_pair(r2, id));
    }
    unsigned n = k-smin.size();
    std::partial_sort(pts.begin(), pts.begin()+n, pts.end());

    /* append these to smin and return */
    for (unsigned i=0; i<n; i++) smin.push_back(pts[i].second);
    std::sort(smin.begin(), smin.end());
    return smin;
}

// do only bounding box checks on query atoms, no spatial hashing.
IdList SpatialHash::find_within_small(float r, const float* pos,
    int n, const Id* ids) const {

    IdList result;
    const float r2 = r*r;
    __m128 rad2 = _mm_set1_ps(r2);
    float xlo = xmin - r;
    float ylo = ymin - r;
    float zlo = zmin - r;

    float xhi = xmax + r;
    float yhi = ymax + r;
    float zhi = zmax + r;

    __m128 lo = _mm_set_ps(0,zlo,ylo,xlo);
    __m128 hi = _mm_set_ps(0,zhi,yhi,xhi);

    for (int i=0; i<n; i++) {
        Id id = ids[i];
        const float* p = pos + 3*id;
        const float x = p[0];
        const float y = p[1];
        const float z = p[2];

        __m128 xyz = _mm_set_ps(0,z,y,x);
        if (_mm_movemask_ps(_mm_cmplt_ps(xyz, lo)) |
            _mm_movemask_ps(_mm_cmpgt_ps(xyz, hi))) {
            continue;
        }

        int j=0;
        __m128 xj = _mm_set1_ps(x);
        __m128 yj = _mm_set1_ps(y);
        __m128 zj = _mm_set1_ps(z);

        for (; j<ntarget-4; j+=4) {
            __m128 p0, p1, p2;
            __m128 q0, q1, q2;

            p0 = _mm_load_ps(_x + j);
            p1 = _mm_load_ps(_y + j);
            p2 = _mm_load_ps(_z + j);

            q0 = _mm_sub_ps(p0, xj);    /* dx */
            q1 = _mm_sub_ps(p1, yj);    /* dy */
            q2 = _mm_sub_ps(p2, zj);    /* dz */
            
            p0 = _mm_mul_ps(q0, q0);    /* dx**2 */
            p1 = _mm_mul_ps(q1, q1);    /* dy**2 */
            p2 = _mm_mul_ps(q2, q2);    /* dz**2 */

            q0 = _mm_add_ps(_mm_add_ps(p0,p1),p2);
            if (_mm_movemask_ps(_mm_cmple_ps(q0, rad2))) {
                result.push_back(id);
                j = ntarget;    // skip the next loop
                break;
            }
        }
        for (; j<ntarget; j++) {
            float dx = x - _x[j];
            float dy = y - _y[j];
            float dz = z - _z[j];
            float d2 = dx*dx + dy*dy + dz*dz;
            if (d2 <= r2) {
                result.push_back(id);
                break;
            }
        }
    }
    return result;
}

IdList SpatialHash::find_within(float r, const float* pos, 
                  int n, const Id* ids) const {

    bool periodic = cx!=0 || cy!=0 || cz!=0;

    // pbwithin not implemented for find_within_small.
    // threshold based on minimal testing.
    if (!periodic && ntarget < 800) return find_within_small(r,pos,n,ids);

    IdList result;
    result.reserve(n);
    int j=0;

    float tmp[12];

#ifdef __SSE2__
   const float r2 = r*r;
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

       const float *a = pos+3*ids[j  ];
       const float *b = pos+3*ids[j+1];
       const float *c = pos+3*ids[j+2];
       const float *d = pos+3*ids[j+3];
       if (rot) {
           memcpy(tmp+0,a,3*sizeof(float));
           memcpy(tmp+3,b,3*sizeof(float));
           memcpy(tmp+6,c,3*sizeof(float));
           memcpy(tmp+9,d,3*sizeof(float));
           pfx::apply_rotation(4,tmp,rot);
           a=tmp+0;
           b=tmp+3;
           c=tmp+6;
           d=tmp+9;
       }
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
        if (periodic) {
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
        const float *xyz = pos + 3*id;
        if (rot) {
            memcpy(tmp,xyz,3*sizeof(float));
            pfx::apply_rotation(1,tmp,rot);
            xyz = tmp;
        }
        float x=xyz[0], y=xyz[1], z=xyz[2];
        if (test(r, x,y,z) || (periodic && minimage(r,cx,cy,cz,x,y,z))) {
            result.push_back(id);
        }
    }
    return result;
}

SpatialHash::ContactList
SpatialHash::findContacts(float r, const float* pos,
                               int n, const Id* ids) {
    //double t=-now();
    contact_array_t arr;
    findContacts(r, pos, n, ids, &arr);
    //t+=now();
    //printf("raw: %.3fms\n", t*1000);
    ContactList result;
    result.reserve(arr.count);
    for (uint64_t i=0, n=arr.count; i<n; i++) {
        result.emplace_back(arr.i[i], arr.j[i], arr.d2[i]);
    }
    return result;
}

void SpatialHash::findContacts(float r, const float* pos,
                               int n, const Id* ids,
                               contact_array_t* result) {
    voxelize(r);
    findContactsReuseVoxels(r, pos, n, ids, result);
}

void SpatialHash::findContactsReuseVoxels(float r, const float* pos,
                                          int n, const Id* ids,
                                          contact_array_t* result) const {

    bool periodic = cx!=0 || cy!=0 || cz!=0;
    float tmp[3];

    for (int j=0; j<n; j++) {
        unsigned id = ids ? ids[j] : j;
        const float *xyz = pos + 3*id;
        if (rot) {
            memcpy(tmp,xyz,3*sizeof(float));
            pfx::apply_rotation(1,tmp,rot);
            xyz = tmp;
        }
        float x=xyz[0], y=xyz[1], z=xyz[2];
        int xi = (x-ox) * ir;
        int yi = (y-oy) * ir;
        int zi = (z-oz) * ir;
        if (!(xi<0 || xi>=nx ||
              yi<0 || yi>=ny ||
              zi<0 || zi>=nz)) {
            int voxid = zi + nz*(yi + ny*xi);
            find_contacts(r*r, voxid, x,y,z, id, result);
        }
        if (periodic) {
            minimage_contacts(r, cx,cy,cz, x,y,z, id, result);
        }
    }
}

void SpatialHash::minimage_contacts(float r, float ga, float gb, float gc,
                                    float px, float py, float pz,
                                    Id id, contact_array_t* result) const {
    float xlo = xmin - r;
    float ylo = ymin - r;
    float zlo = zmin - r;
    float xhi = xmax + r;
    float yhi = ymax + r;
    float zhi = zmax + r;
    for (int i=-1; i<=1; i++) {
        float x = px + ga*i;
        if (x<xlo || x>xhi) continue;
        for (int j=-1; j<=1; j++) {
            float y = py + gb*j;
            if (y<ylo || y>yhi) continue;
            for (int k=-1; k<=1; k++) {
                float z = pz + gc*k;
                if (z<zlo || z>zhi) continue;
                if (i==0 && j==0 && k==0) continue;
                int xi = (x-ox) * ir;
                int yi = (y-oy) * ir;
                int zi = (z-oz) * ir;
                if (xi<0 || xi>=nx ||
                    yi<0 || yi>=ny ||
                    zi<0 || zi>=nz) continue;
                int voxid = zi + nz*(yi + ny*xi);
                find_contacts(r*r, voxid, x,y,z, id, result);
            }
        }
    }
}

bool SpatialHash::minimage(float r, float ga, float gb, float gc,
                           float px, float py, float pz) const {
    float xlo = xmin - r;
    float ylo = ymin - r;
    float zlo = zmin - r;
    float xhi = xmax + r;
    float yhi = ymax + r;
    float zhi = zmax + r;
    for (int i=-1; i<=1; i++) {
        float x = px + ga*i;
        if (x<xlo || x>xhi) continue;
        for (int j=-1; j<=1; j++) {
            float y = py + gb*j;
            if (y<ylo || y>yhi) continue;
            for (int k=-1; k<=1; k++) {
                float z = pz + gc*k;
                if (z<zlo || z>zhi) continue;
                if (i==0 && j==0 && k==0) continue;
                if (test(r,x,y,z)) return true;
            }
        }
    }
    return false;
}

static const uint64_t lut_zero = 0x8080808080808080;
static const uint64_t lutvals[] = {
    lut_zero,              lut_zero,  // 0
    0x8080808003020100,    lut_zero,  // 1
    0x8080808007060504,    lut_zero,  // 2
    0x0706050403020100,    lut_zero,  // 3
    0x808080800b0a0908,    lut_zero,  // 4
    0x0b0a090803020100,    lut_zero,  // 5
    0x0b0a090807060504,    lut_zero,  // 6    0110
    0x0706050403020100,    0x808080800b0a0908,  // 7    0111
    0x808080800f0e0d0c,    lut_zero,  // 8
    0x0f0e0d0c03020100,    lut_zero,  // 9    1001
    0x0f0e0d0c07060504,    lut_zero,  // 10   1010
    0x0706050403020100,    0x808080800f0e0d0c,  // 11   1011
    0x0f0e0d0c0b0a0908,    lut_zero,  // 12   1100
    0x0b0a090803020100,    0x808080800f0e0d0c,  // 13   1101
    0x0b0a090807060504,    0x808080800f0e0d0c,  // 14   1110
    0x0706050403020100,    0x0f0e0d0c0b0a0908,  // 15   1111
};
static const __m128i* lut = (const __m128i*)(lutvals);

// popcnt for 0-15
static const uint8_t popcnt_u4_data[] = {
    0,1,1,2,
    1,2,2,3,
    1,2,2,3,
    2,3,3,4
};
static inline uint8_t my_popcnt_u4(int mask) {
    return popcnt_u4_data[mask];
}

void SpatialHash::find_contacts(float r2, int voxid, float x, float y, float z,
                           Id id, contact_array_t* result) const {

    result->reserve_additional(27*maxcount);
    Id* ri = result->i;
    Id* rj = result->j;
    float* rd = result->d2;
    uint64_t count = result->count;

#ifdef __SSE4_1__
    __m128i packed_i = _mm_set1_epi32(id);
    __m128 xj = _mm_set1_ps(x);
    __m128 yj = _mm_set1_ps(y);
    __m128 zj = _mm_set1_ps(z);
    __m128 R2 = _mm_set1_ps(r2);
#endif
    for (int i=0; i<10; i++) {
        int vox = voxid + full_shell[i];
        uint32_t b = _counts[vox], e = _counts[vox+strip_lens[i]];
        const float * xi = _x+b;
        const float * yi = _y+b;
        const float * zi = _z+b;
        const Id    * ii = _ids+b;

#ifdef __SSE4_1__
        /* advance to aligned offset */
        for (; b<e && (b&3); ++b, ++xi, ++yi, ++zi, ++ii) {
#ifdef EXCLUDE_SELF_CONTACTS
            if (id==*ii) continue;
#endif
            float dx = x - *xi;
            float dy = y - *yi;
            float dz = z - *zi;
            float d2 = dx*dx + dy*dy + dz*dz;
            if (d2<=r2) {
                ri[count] = id;
                rj[count] = *ii;
                rd[count] = d2;
                ++count;
            }
        }

        /* simd for chunks of four */
        for (; b+4<e; b+=4, xi+=4, yi+=4, zi+=4, ii+=4) {
            __m128 p0, p1, p2;
            __m128 q0, q1, q2;
            __m128i i4;

            p0 = _mm_load_ps(xi);
            p1 = _mm_load_ps(yi);
            p2 = _mm_load_ps(zi);
            i4 = _mm_load_si128((const __m128i*)ii);

            q0 = _mm_sub_ps(p0, xj);    /* dx */
            q1 = _mm_sub_ps(p1, yj);    /* dy */
            q2 = _mm_sub_ps(p2, zj);    /* dz */
            
            p0 = _mm_mul_ps(q0, q0);    /* dx**2 */
            p1 = _mm_mul_ps(q1, q1);    /* dy**2 */
            p2 = _mm_mul_ps(q2, q2);    /* dz**2 */

            q0 = _mm_add_ps(_mm_add_ps(p0,p1),p2);  // d2

            // mask contains a value in the range [0,15] whose bits 
            // correspond to the values in range.
            int mask = _mm_movemask_ps(
#ifdef EXCLUDE_SELF_CONTACTS
                    _mm_andnot_ps(
                        _mm_castsi128_ps(_mm_cmpeq_epi32(i4, packed_i)),
                        _mm_cmple_ps(q0, R2)));
#else
                    _mm_cmple_ps(q0, R2));
#endif


            // push the kept indices and distances into a contiguous chunk
            __m128i shuf = lut[mask];
            __m128i packed_d2 = _mm_shuffle_epi8(_mm_castps_si128(q0), shuf);
            __m128i packed_j = _mm_shuffle_epi8(i4, shuf);

            // write to destination arrays and update count
            _mm_storeu_si128((__m128i*)(ri+count), packed_i);
            _mm_storeu_si128((__m128i*)(rj+count), packed_j);
            _mm_storeu_si128((__m128i*)(rd+count), packed_d2);
            //count += _mm_popcnt_u32(mask);    // needs sse4_2
            count += my_popcnt_u4(mask);
        }
#endif
        /* stragglers */
        for (; b<e; ++b, ++xi, ++yi, ++zi, ++ii) {
#ifdef EXCLUDE_SELF_CONTACTS
            if (id==*ii) continue;
#endif
            float dx = x - *xi;
            float dy = y - *yi;
            float dz = z - *zi;
            float d2 = dx*dx + dy*dy + dz*dz;
            if (d2<=r2) {
                ri[count] = id;
                rj[count] = *ii;
                rd[count] = d2;
                ++count;
            }
        }
    }
    result->count = count;
}


bool SpatialHash::test2(float r2, int voxid, float x, float y, float z) const {
#ifdef __SSE2__
    __m128 xj = _mm_set1_ps(x);
    __m128 yj = _mm_set1_ps(y);
    __m128 zj = _mm_set1_ps(z);
    __m128 R2 = _mm_set1_ps(r2);
#endif
    for (int i=0; i<10; i++) {
        int vox = voxid + full_shell[i];
        uint32_t b = _counts[vox], e = _counts[vox+strip_lens[i]];
        const float * xi = _x+b;
        const float * yi = _y+b;
        const float * zi = _z+b;

#ifdef __SSE2__
        /* advance to aligned offset */
        for (; b<e && (b&3); ++b, ++xi, ++yi, ++zi) {
            float dx = x - *xi;
            float dy = y - *yi;
            float dz = z - *zi;
            float d2 = dx*dx + dy*dy + dz*dz;
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
            float dx = x - *xi;
            float dy = y - *yi;
            float dz = z - *zi;
            float d2 = dx*dx + dy*dy + dz*dz;
            if (d2<=r2) return true;
        }
    }
    return false;
}


namespace desres { namespace msys {
    IdList FindWithin(IdList const& wsel, const float* wat,
                      IdList const& psel, const float* pro,
                      float radius,
                      const double* cell) {
 
        return SpatialHash(pro, psel.size(), &psel[0], cell)
            .findWithin(radius, wat, wsel.size(), &wsel[0]);
    }
 
    IdList FindNearest(IdList const& wsel, const float* wat,
                       IdList const& psel, const float* pro,
                       Id k,               const double* cell) {
 
        return SpatialHash(pro, psel.size(), &psel[0], cell)
            .findNearest(k, wat, wsel.size(), &wsel[0]);
    }

}}


