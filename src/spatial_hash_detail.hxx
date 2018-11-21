#include <assert.h>

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

template <typename Float>
void SpatialHashT<Float>::compute_full_shell() {
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

template <typename Float>
SpatialHashT<Float>::~SpatialHashT() {
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

template <typename Float>
SpatialHashT<Float>::SpatialHashT( const Float *pos, int n, const Id* ids, 
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
#ifdef WIN32
        rot = (Float*)_aligned_malloc(9*sizeof(*rot), 16);
#else
        assert(0==posix_memalign((void **)&rot, 16, 9*sizeof(*rot)));
#endif
        double d1=0, d2=0, d3=0; /* row dot-products */
        for (int i=0; i<3; i++) {
            rot[0+i] = cell[0+i]/cx;
            rot[3+i] = cell[3+i]/cy;
            rot[6+i] = cell[6+i]/cz;
            d1 += rot[0+i]*rot[3+i];
            d2 += rot[0+i]*rot[6+i];
            d3 += rot[3+i]*rot[6+i];
        }
        static const Float eps = 1e-4;
        if (fabs(d1)>eps || fabs(d2)>eps || fabs(d3)>eps) {
            MSYS_FAIL("cell appears triclinic: dot products " << d1 << " " << d2 << " " << d3);
        }
        if (!(rot[1] || rot[2] || rot[3] || rot[5] || rot[6] || rot[7])) {
            free(rot);
            rot=NULL;
        }
    }

    /* copy to transposed arrays */
#ifdef WIN32
    _x = (Float*)_aligned_malloc(ntarget*sizeof(*_x), 16);
    _y = (Float*)_aligned_malloc(ntarget*sizeof(*_y), 16);
    _z = (Float*)_aligned_malloc(ntarget*sizeof(*_z), 16);
    _tmpx = (Float*)_aligned_malloc(ntarget*sizeof(*_tmpx), 16);
    _tmpy = (Float*)_aligned_malloc(ntarget*sizeof(*_tmpy), 16);
    _tmpz = (Float*)_aligned_malloc(ntarget*sizeof(*_tmpz), 16);
#else
    assert(0==posix_memalign((void **)&_x, 16, ntarget*sizeof(*_x)));
    assert(0==posix_memalign((void **)&_y, 16, ntarget*sizeof(*_y)));
    assert(0==posix_memalign((void **)&_z, 16, ntarget*sizeof(*_z)));
    assert(0==posix_memalign((void **)&_tmpx, 16, ntarget*sizeof(*_tmpx)));
    assert(0==posix_memalign((void **)&_tmpy, 16, ntarget*sizeof(*_tmpy)));
    assert(0==posix_memalign((void **)&_tmpz, 16, ntarget*sizeof(*_tmpz)));
    assert(0==posix_memalign((void **)&_ids, 16, ntarget*sizeof(*_ids)));
    assert(0==posix_memalign((void **)&_tmpids, 16, ntarget*sizeof(*_tmpids)));
#endif
    for (int i=0; i<ntarget; i++) {
        Id id = ids ? ids[i] : i;
        const Float *xyz = pos+3*id;
        if(rot) {
            Float pos[3]={xyz[0],xyz[1],xyz[2]};
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


template <typename Float>
SpatialHashT<Float>& SpatialHashT<Float>::voxelize(Float r) {
    if (r<=0) MSYS_FAIL("radius " << r << " must be positive");
    rad = r;
    ir = Float(1)/rad;
    ox = xmin - rad;
    oy = ymin - rad;
    oz = zmin - rad;

    /* construct voxel grid.  */
    nx = (xmax-xmin)*ir + 3;
    ny = (ymax-ymin)*ir + 3;
    nz = (zmax-zmin)*ir + 3;
    static const int maxdim = 500;
    if (nx > maxdim || ny > maxdim || nz > maxdim) {
        Float dbig = std::max(std::max(xmax-xmin, ymax-ymin), zmax-zmin);
        ir = Float(maxdim) / dbig;
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
        Float x = _x[i];
        Float y = _y[i];
        Float z = _z[i];
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

template <typename Float>
Float SpatialHashT<Float>::mindist2(Float x, Float y, Float z) const {
    int xi = (x-ox) * ir;
    int yi = (y-oy) * ir;
    int zi = (z-oz) * ir;
    Float r2 = std::numeric_limits<Float>::max();
    if (xi<0 || xi>=nx ||
        yi<0 || yi>=ny ||
        zi<0 || zi>=nz) return r2;
    int voxid = zi + nz*(yi + ny*xi);
    /* TODO: SIMD */
    for (int i=0; i<10; i++) {
        int vox = voxid + full_shell[i];
        uint32_t b = _counts[vox], e = _counts[vox+strip_lens[i]];
        const Float *xi = _x+b;
        const Float *yi = _y+b;
        const Float *zi = _z+b;

        for (; b<e; ++b, ++xi, ++yi, ++zi) {
            Float dx = x - *xi;
            Float dy = y - *yi;
            Float dz = z - *zi;
            Float d2 = dx*dx + dy*dy + dz*dz;
            r2 = std::min(r2, d2);
        }
    }
    return r2;
}

template <typename Float>
IdList SpatialHashT<Float>::findNearest(unsigned k, const Float* pos, 
                   unsigned num, const Id* ids) {

    /* If there are no hashed atoms, we should stop now. */
    if (ntarget==0) MSYS_FAIL("No atoms in target selection");

    /* enough atoms available */
    if (num <= k) return IdList(ids, ids+num);

    IdList smin, smax;
    Float rmin = 0.0;
    Float rmax = 2.5;

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
        Float rm = 0.5*(rmin+rmax);
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
    std::vector<std::pair<Float, Id> > pts;
    pts.reserve(smax.size()-smin.size());
    for (unsigned i=0, n=smax.size(); i<n; i++) {
        Id id = smax[i];
        if (std::binary_search(smin.begin(), smin.end(), id)) continue;
        const Float* p = pos+3*id;
        Float r2 = mindist2(p[0], p[1], p[2]);
        pts.push_back(std::make_pair(r2, id));
    }
    unsigned n = k-smin.size();
    std::partial_sort(pts.begin(), pts.begin()+n, pts.end());

    /* append these to smin and return */
    for (unsigned i=0; i<n; i++) smin.push_back(pts[i].second);
    std::sort(smin.begin(), smin.end());
    return smin;
}

template <typename Float>
typename SpatialHashT<Float>::ContactList
SpatialHashT<Float>::findContacts(Float r, const Float* pos,
                               int n, const Id* ids) {
    contact_array_t arr;
    findContacts(r, pos, n, ids, &arr);
    ContactList result;
    result.reserve(arr.count);
    for (uint64_t i=0, n=arr.count; i<n; i++) {
        result.emplace_back(arr.i[i], arr.j[i], arr.d2[i]);
    }
    return result;
}

template <typename Float>
void SpatialHashT<Float>::findContacts(Float r, const Float* pos,
                               int n, const Id* ids,
                               contact_array_t* result) {
    voxelize(r);
    findContactsReuseVoxels(r, pos, n, ids, result);
}

template <typename Float>
void SpatialHashT<Float>::findContactsReuseVoxels(Float r, const Float* pos,
                                          int n, const Id* ids,
                                          contact_array_t* result) const {

    bool periodic = cx!=0 || cy!=0 || cz!=0;
    Float tmp[3];

    for (int j=0; j<n; j++) {
        unsigned id = ids ? ids[j] : j;
        const Float *xyz = pos + 3*id;
        if (rot) {
            memcpy(tmp,xyz,3*sizeof(Float));
            pfx::apply_rotation(1,tmp,rot);
            xyz = tmp;
        }
        Float x=xyz[0], y=xyz[1], z=xyz[2];
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
template <typename Float>
template <typename SpatialHashExclusions>
void SpatialHashT<Float>::findPairlistReuseVoxels(Float r, SpatialHashExclusions const& excl, contact_array_t *result) const {
    bool periodic = cx!=0 || cy!=0 || cz!=0;
    Float tmp[3];

    for (int j=0; j<ntarget; j++) {
        unsigned id = _ids[j];
        Float x = _x[j];
        Float y = _y[j];
        Float z = _z[j];
        if (rot) {
            tmp[0] = x;
            tmp[1] = y;
            tmp[2] = z;
            pfx::apply_rotation(1,tmp,rot);
            x = tmp[0];
            y = tmp[1];
            z = tmp[2];
        }
        int xi = (x-ox) * ir;
        int yi = (y-oy) * ir;
        int zi = (z-oz) * ir;
        if (!(xi<0 || xi>=nx ||
              yi<0 || yi>=ny ||
              zi<0 || zi>=nz)) {
            int voxid = zi + nz*(yi + ny*xi);
            find_pairlist(r*r, voxid, x,y,z, id, excl, result);
        }
        if (periodic) {
            minimage_pairlist(r, cx,cy,cz, x,y,z, id, excl, result);
        }
    }
}

template <typename Float>
template <typename SpatialHashExclusions>
void SpatialHashT<Float>::find_pairlist(Float r2, int voxid, Float x, Float y, Float z,
                           Id id, SpatialHashExclusions const& excl, contact_array_t* result) const {

    result->reserve_additional(27*maxcount);
    Id* ri = result->i;
    Id* rj = result->j;
    Float* rd = result->d2;
    uint64_t count = result->count;
    uint64_t key_hi(id);
    key_hi <<= 32;

    for (int i=0; i<10; i++) {
        int vox = voxid + full_shell[i];
        uint32_t b = _counts[vox], e = _counts[vox+strip_lens[i]];
        const Float * xi = _x+b;
        const Float * yi = _y+b;
        const Float * zi = _z+b;
        const Id    * ii = _ids+b;

        /* stragglers */
        for (; b<e; ++b, ++xi, ++yi, ++zi, ++ii) {
          if (id >= *ii) continue;
            Float dx = x - *xi;
            Float dy = y - *yi;
            Float dz = z - *zi;
            Float d2 = dx*dx + dy*dy + dz*dz;
            if (d2<=r2 && !excl.count(key_hi | *ii)) {
                ri[count] = id;
                rj[count] = *ii;
                rd[count] = d2;
                ++count;
            }
        }
    }
    result->count = count;
}

template <typename Float>
void SpatialHashT<Float>::minimage_contacts(Float r, Float ga, Float gb, Float gc,
                                    Float px, Float py, Float pz,
                                    Id id, contact_array_t* result) const {
    Float xlo = xmin - r;
    Float ylo = ymin - r;
    Float zlo = zmin - r;
    Float xhi = xmax + r;
    Float yhi = ymax + r;
    Float zhi = zmax + r;
    for (int i=-1; i<=1; i++) {
        Float x = px + ga*i;
        if (x<xlo || x>xhi) continue;
        for (int j=-1; j<=1; j++) {
            Float y = py + gb*j;
            if (y<ylo || y>yhi) continue;
            for (int k=-1; k<=1; k++) {
                Float z = pz + gc*k;
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

template <typename Float>
template <typename SpatialHashExclusions>
void SpatialHashT<Float>::minimage_pairlist(Float r, Float ga, Float gb, Float gc,
                                            Float px, Float py, Float pz,
                                            Id id, SpatialHashExclusions const& excl,
                                            contact_array_t* result) const {
    Float xlo = xmin - r;
    Float ylo = ymin - r;
    Float zlo = zmin - r;
    Float xhi = xmax + r;
    Float yhi = ymax + r;
    Float zhi = zmax + r;
    for (int i=-1; i<=1; i++) {
        Float x = px + ga*i;
        if (x<xlo || x>xhi) continue;
        for (int j=-1; j<=1; j++) {
            Float y = py + gb*j;
            if (y<ylo || y>yhi) continue;
            for (int k=-1; k<=1; k++) {
                Float z = pz + gc*k;
                if (z<zlo || z>zhi) continue;
                if (i==0 && j==0 && k==0) continue;
                int xi = (x-ox) * ir;
                int yi = (y-oy) * ir;
                int zi = (z-oz) * ir;
                if (xi<0 || xi>=nx ||
                    yi<0 || yi>=ny ||
                    zi<0 || zi>=nz) continue;
                int voxid = zi + nz*(yi + ny*xi);
                find_pairlist(r*r, voxid, x,y,z, id, excl, result);
            }
        }
    }
}

template <typename Float>
bool SpatialHashT<Float>::minimage(Float r, Float ga, Float gb, Float gc,
                           Float px, Float py, Float pz) const {
    Float xlo = xmin - r;
    Float ylo = ymin - r;
    Float zlo = zmin - r;
    Float xhi = xmax + r;
    Float yhi = ymax + r;
    Float zhi = zmax + r;
    for (int i=-1; i<=1; i++) {
        Float x = px + ga*i;
        if (x<xlo || x>xhi) continue;
        for (int j=-1; j<=1; j++) {
            Float y = py + gb*j;
            if (y<ylo || y>yhi) continue;
            for (int k=-1; k<=1; k++) {
                Float z = pz + gc*k;
                if (z<zlo || z>zhi) continue;
                if (i==0 && j==0 && k==0) continue;
                if (test(r,x,y,z)) return true;
            }
        }
    }
    return false;
}


