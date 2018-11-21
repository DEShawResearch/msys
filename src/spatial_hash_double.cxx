#include "spatial_hash.hxx"
#include <limits>
#include <stdio.h>

#define EXCLUDE_SELF_CONTACTS

namespace desres { namespace msys {

// do only bounding box checks on query atoms, no spatial hashing.
template<>
IdList SpatialHashT<double>::find_within_small(double r, const double* pos,
    int n, const Id* ids) const {

    IdList result;
    const double r2 = r*r;
    double xlo = xmin - r;
    double ylo = ymin - r;
    double zlo = zmin - r;

    double xhi = xmax + r;
    double yhi = ymax + r;
    double zhi = zmax + r;

    for (int i=0; i<n; i++) {
        Id id = ids[i];
        const double* p = pos + 3*id;
        const double x = p[0];
        const double y = p[1];
        const double z = p[2];

        if (x<xlo || y<ylo || z<zlo ||
            x>xhi || y>yhi || z>zhi) {
            continue;
        }

        for (int j=0; j<ntarget; j++) {
            double dx = x - _x[j];
            double dy = y - _y[j];
            double dz = z - _z[j];
            double d2 = dx*dx + dy*dy + dz*dz;
            if (d2 <= r2) {
                result.push_back(id);
                break;
            }
        }
    }
    return result;
}

template<>
bool SpatialHashT<double>::test2(double r2, int voxid, double x, double y, double z) const {
    for (int i=0; i<10; i++) {
        int vox = voxid + full_shell[i];
        uint32_t b = _counts[vox], e = _counts[vox+strip_lens[i]];
        const double * xi = _x+b;
        const double * yi = _y+b;
        const double * zi = _z+b;

        for (; b<e; ++b, ++xi, ++yi, ++zi) {
            double dx = x - *xi;
            double dy = y - *yi;
            double dz = z - *zi;
            double d2 = dx*dx + dy*dy + dz*dz;
            if (d2<=r2) return true;
        }
    }
    return false;
}

template<>
IdList SpatialHashT<double>::find_within(double r, const double* pos, 
                  int n, const Id* ids) const {

    bool periodic = cx!=0 || cy!=0 || cz!=0;

    // pbwithin not implemented for find_within_small.
    // threshold based on minimal testing.
    if (!periodic && ntarget < 800) return find_within_small(r,pos,n,ids);

    IdList result;
    result.reserve(n);
    int j=0;

    double tmp[12];

    for (; j<n; j++) {
        unsigned id = ids[j];
        const double *xyz = pos + 3*id;
        if (rot) {
            memcpy(tmp,xyz,3*sizeof(double));
            pfx::apply_rotation(1,tmp,rot);
            xyz = tmp;
        }
        double x=xyz[0], y=xyz[1], z=xyz[2];
        if (test(r, x,y,z) || (periodic && minimage(r,cx,cy,cz,x,y,z))) {
            result.push_back(id);
        }
    }
    return result;
}

template<>
void SpatialHashT<double>::find_contacts(double r2, int voxid, double x, double y, double z,
                           Id id, contact_array_t* result) const {

    result->reserve_additional(27*maxcount);
    Id* ri = result->i;
    Id* rj = result->j;
    double* rd = result->d2;
    uint64_t count = result->count;

    for (int i=0; i<10; i++) {
        int vox = voxid + full_shell[i];
        uint32_t b = _counts[vox], e = _counts[vox+strip_lens[i]];
        const double * xi = _x+b;
        const double * yi = _y+b;
        const double * zi = _z+b;
        const Id    * ii = _ids+b;

        /* stragglers */
        for (; b<e; ++b, ++xi, ++yi, ++zi, ++ii) {
#ifdef EXCLUDE_SELF_CONTACTS
            if (id==*ii) continue;
#endif
            double dx = x - *xi;
            double dy = y - *yi;
            double dz = z - *zi;
            double d2 = dx*dx + dy*dy + dz*dz;
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

}}

