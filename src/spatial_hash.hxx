#ifndef desres_msys_spatial_hash_hxx
#define desres_msys_spatial_hash_hxx

#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <vector>
#include <limits>
#include <algorithm>
#include "types.hxx"

namespace desres { namespace msys {

    class SpatialHash {

        /* Support distance queries of at most _radius */
        float rad;
        
        /* inverse radius */
        float ir;

        /* bounding box of target coordinates */
        float xmin, ymin, zmin;
        float xmax, ymax, zmax;

        /* location of grid origin */
        float ox, oy, oz;

        /* voxel grid dimension */
        int nx, ny, nz;

        /* neighbor full shell offsets */
        int full_shell[10];
        int strip_lens[10];

        /* target coordinates */
        const int ntarget;
        float *_x, *_y, *_z;
        float *_tmpx, *_tmpy, *_tmpz;

        /* cumulative number of particles hashed to the i'th voxel */
        std::vector<uint32_t> _countspace;
        uint32_t *_counts;

        void compute_full_shell();

    public:
        ~SpatialHash();
        SpatialHash( const float *pos, int n, const Id* ids);

        SpatialHash& voxelize(float r);

        /* Is the given point within r of any point in the region? */
        bool test(float r, float x, float y, float z) const {
            int xi = (x-ox) * ir;
            int yi = (y-oy) * ir;
            int zi = (z-oz) * ir;
            if (xi<0 || xi>=nx ||
                yi<0 || yi>=ny ||
                zi<0 || zi>=nz) return false;
            int voxid = zi + nz*(yi + ny*xi);
            return test2(r*r, voxid, x,y,z);
        }

        float mindist2(float x, float y, float z) const;

        IdList findNearest(unsigned k, const float* pos, 
                           unsigned num, const Id* ids,
                           const double* cell = NULL);

        IdList findWithin(float r, const float* pos,
                          int n, const Id* ids,
                          const double* cell = NULL) {
            voxelize(r);
            return find_within(r,pos,n,ids,cell);
        }

        IdList find_within(float r, const float* pos, 
                          int n, const Id* ids,
                          const double* cell ) const;

        bool minimage(float r, float ga, float gb, float gc,
                                float px, float py, float pz) const;
        bool test2(float r2, int voxid, float x, float y, float z) const;
        float radius() const { return rad; }
    };

    /* convenience routines */
    IdList FindWithin(IdList const& wsel, const float* wat,
                      IdList const& psel, const float* pro,
                      float radius,
                      const double* cell);
 
    IdList FindNearest(IdList const& wsel, const float* wat,
                       IdList const& psel, const float* pro,
                       Id k,               const double* cell);
 
}}

#endif

