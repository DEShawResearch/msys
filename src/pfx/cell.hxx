#ifndef desres_pfx_cell_hxx
#define desres_pfx_cell_hxx

namespace desres { namespace msys { namespace pfx {

    /* Construct new cell with each vector scaled by 1/|v|^2. */
    template <typename scalar>
    void make_projection(const scalar* cell, scalar* proj) {
        for (int i=0; i<3; i++) {
            scalar x,y,z,n2;
            x = *cell++;
            y = *cell++;
            z = *cell++;
            n2 = x*x + y*y + z*z;
            if (n2>0) {
                n2 = ((scalar)1)/n2;
            }
            *proj++ = x*n2;
            *proj++ = y*n2;
            *proj++ = z*n2;
        }
    }

    inline double ROUND(double x) {
        static const double BLACK_MAGIC = 6755399441055744.0;
        // need to write out to memory to force loss of precision
        volatile double temp = x + BLACK_MAGIC;
        return temp - BLACK_MAGIC;
    }

    inline float ROUND(float x) {
        static const float BLACK_MAGIC = 12582912.0;
        // need to write out to memory to force loss of precision
        volatile float temp = x + BLACK_MAGIC;
        return temp - BLACK_MAGIC;
    }

    template <typename scalar>
    void wrap_vector(const scalar* cell, const scalar* proj,
                     scalar* pos) {
        scalar nx, ny, nz;
        scalar x = pos[0];
        scalar y = pos[1];
        scalar z = pos[2];
    
        nx = -ROUND(proj[0]*x + proj[1]*y + proj[2]*z);
        ny = -ROUND(proj[3]*x + proj[4]*y + proj[5]*z);
        nz = -ROUND(proj[6]*x + proj[7]*y + proj[8]*z);
    
        x = nx*cell[0] + ny*cell[3] + nz*cell[6];
        y = nx*cell[1] + ny*cell[4] + nz*cell[7];
        z = nx*cell[2] + ny*cell[5] + nz*cell[8];
    
        pos[0]=x;
        pos[1]=y;
        pos[2]=z;
    }

    template <typename scalar, typename center_scalar, typename wscalar>
    void compute_center( unsigned nelems,
                         const unsigned* elems,
                         const scalar* pos,
                         center_scalar *center,
                         const wscalar* wts = nullptr) {
        scalar cx=0, cy=0, cz=0, cw=0;
        for (unsigned i=0; i<nelems; i++) {
            const scalar *p = pos + (elems ? 3*elems[i] : 3*i);
            const scalar w = wts ? *wts++ : 1;
            cw += w;
            cx += w*p[0];
            cy += w*p[1];
            cz += w*p[2];
        }
        cw = cw ? 1/cw : 0;
        center[0] = cx * cw;
        center[1] = cy * cw;
        center[2] = cz * cw;
    }

    template <typename scalar>
    void apply_shift(unsigned n,
                     scalar* pos,
                     scalar cx,
                     scalar cy,
                     scalar cz) {
        for (; n>0; n--) {
            *pos++ += cx;
            *pos++ += cy;
            *pos++ += cz;
        }
    }


}}}

#endif
