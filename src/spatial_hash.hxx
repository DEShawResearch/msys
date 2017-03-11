#ifndef desres_msys_spatial_hash_hxx
#define desres_msys_spatial_hash_hxx

#include "types.hxx"
#include <math.h>
#include "pfx/rms.hxx"

namespace desres { namespace msys {

    template <typename Float>
    class SpatialHashT {

        /* Support distance queries of at most _radius */
        Float rad;
        
        /* inverse voxel size */
        Float ir;

        /* bounding box of target coordinates */
        Float xmin, ymin, zmin;
        Float xmax, ymax, zmax;

        /* location of grid origin */
        Float ox, oy, oz;

        /* voxel grid dimension */
        int nx, ny, nz;

        /* rotation */
        Float* rot;
        /* cell dims */
        Float cx, cy, cz;

        /* neighbor full shell offsets */
        int full_shell[10];
        int strip_lens[10];

        /* target coordinates */
        const int ntarget;
        Float *_x, *_y, *_z;
        Float *_tmpx, *_tmpy, *_tmpz;

        /* cumulative number of particles hashed to the i'th voxel */
        std::vector<uint32_t> _countspace;
        uint32_t *_counts;
        uint32_t maxcount=0;

        /* ids of hashed particles */
        Id* _ids;
        Id *_tmpids;

        void compute_full_shell();
        bool test2(Float r2, int voxid, Float x, Float y, Float z) const;

        static
        void find_bbox(int n, const Float* x, Float *_min, Float *_max) {
            Float min = x[0], max=x[0];
            for (int i=0; i<n; i++) {
                Float p = x[i];
                min = min < p ? min : p;
                max = max > p ? max : p;
            }
            *_min = min;
            *_max = max;
        }

    public:
        ~SpatialHashT();

        /* Constructor: supply n ids of points to be hashed. */
        SpatialHashT(const Float *pos, int n, const Id* ids, const double* cell);

        SpatialHashT& voxelize(Float r);

        /* Is the given point within r of any point in the region? */
        bool test(Float r, Float x, Float y, Float z) const {
            int xi = (x-ox) * ir;
            int yi = (y-oy) * ir;
            int zi = (z-oz) * ir;
            if (xi<0 || xi>=nx ||
                yi<0 || yi>=ny ||
                zi<0 || zi>=nz) return false;
            int voxid = zi + nz*(yi + ny*xi);
            return test2(r*r, voxid, x,y,z);
        }

        /* minimum square distance to the the given point from a 
         * hashed point. */
        Float mindist2(Float x, Float y, Float z) const;

        /* find the ids of the k nearest points.  */
        IdList findNearest(unsigned k, const Float* pos, 
                           unsigned num, const Id* ids);

        /* Find the points within r of the hashed points.  If cell
         * is not NULL, use minimum image distances.  */
        IdList findWithin(Float r, const Float* pos,
                          int n, const Id* ids) {
            voxelize(r);
            return find_within(r,pos,n,ids);
        }

        struct contact_t {
            Id i;
            Id j;
            Float d2;
            contact_t(Id i, Id j, Float d2) : i(i), j(j), d2(d2) {}
            bool operator==(contact_t const& c) const {
                return i==c.i && j==c.j && d2==c.d2;
            }
            bool operator<(contact_t const& c) const {
                if (i!=c.i) return i<c.i;
                return j<c.j;
            }
        };
        typedef std::vector<contact_t> ContactList;
        ContactList findContacts(Float r, const Float* pos,
                                 int n, const Id* ids);

        struct contact_array_t {
            Id* i = nullptr;
            Id* j = nullptr;
            Float* d2 = nullptr;
            uint64_t count = 0;
            uint64_t max_size = 0;
            
            ~contact_array_t() {
                free(i);
                free(j);
                free(d2);
            }
            contact_array_t() {}
            contact_array_t(contact_array_t const&) = delete;
            contact_array_t& operator=(contact_array_t const&) = delete;

            void reserve_additional(uint64_t extra) {
                uint64_t size = count + extra;
                if (size > max_size) {
                    max_size = std::max(max_size + max_size/2, size);
                    i = (Id*)realloc(i, max_size*sizeof(*i));
                    j = (Id*)realloc(j, max_size*sizeof(*j));
                    d2= (Float*)realloc(d2,max_size*sizeof(*d2));
                }
            }
        };

        /* find contacts, performing a voxelization step before
         * the search (making this a non-const method)
         */
        void findContacts(Float r, const Float* pos,
                          int n, const Id* ids,
                          contact_array_t* result);

        /* find contacts assuming voxelize(R) for for R>=r has already been
         * called.  const and reentrant.  */
        void findContactsReuseVoxels(Float r, const Float* pos,
                                     int n, const Id* ids,
                                     contact_array_t* result) const;

        /* For expert users only.  Finds points within r of the
         * hashed points assuming the hashed points have already
         * been voxelized with a grid spacing of at least r. */
        IdList find_within(Float r, const Float* pos, 
                          int n, const Id* ids) const;

        IdList find_within_small(Float r, const Float* pos, 
                                 int n, const Id* ids) const;
 
        void find_contacts(Float r2, int voxid, Float x, Float y, Float z,
                           Id id, contact_array_t* result) const;

        void minimage_contacts(Float r, Float ga, Float gb, Float gc,
                               Float px, Float py, Float pz,
                               Id id, contact_array_t* result) const;

        /* Return true if point px,py,pz is within r of some hashed
         * point assuming an orthorhombic periodic cell with lengths
         * ga,gb,gc. */
        bool minimage(Float r, Float ga, Float gb, Float gc,
		      Float px, Float py, Float pz) const;

        /* last voxelization radius */
        Float radius() const { return rad; }
    };

    typedef SpatialHashT<float> SpatialHash;

    /* convenience routines */
    inline IdList FindWithin(IdList const& wsel, const float* wat,
                      IdList const& psel, const float* pro,
                      float radius,
                      const double* cell) {
        return SpatialHash(pro, psel.size(), &psel[0], cell)
            .findWithin(radius, wat, wsel.size(), &wsel[0]);
    }
 
    inline IdList FindNearest(IdList const& wsel, const float* wat,
                       IdList const& psel, const float* pro,
                       Id k,               const double* cell) {
        return SpatialHash(pro, psel.size(), &psel[0], cell)
            .findNearest(k, wat, wsel.size(), &wsel[0]);
    }
 
#include "spatial_hash_detail.hxx"

}}

#endif

