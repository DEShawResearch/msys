#ifndef desres_pfx_rms_hxx
#define desres_pfx_rms_hxx

#include "svd.hxx"
#include <cmath>

namespace desres { namespace msys { namespace pfx {

    /* 3x3 matrix transpose */
    inline void trans_3x3(double * a_tr, const double * a) {
        int i,j;
        for(i=0; i<3; ++i) 
            for(j=0; j<3; ++j) 
                a_tr[3*i + j] = a[3*j + i];
    }
    
    //*3x3 matrix-matrix multiply */
    inline void matmult_3x3(double* c, const double* a, const double* b) {
        int i,j,k;
        for(i=0; i<3; ++i) {
            for(j=0; j<3; ++j) {
                c[3*i + j] = 0;
                for(k=0; k<3; ++k) 
                    c[3*i + j] += a[3*i + k] * b[3*k + j];
            }
        }
    }

    template <typename scalar>
    inline void apply_rotation(int n, scalar* pos, const scalar* mat) {
#define M(i,j) (mat[3*i+j])
        for (;n>0; n--) {
            const scalar x = pos[0];
            const scalar y = pos[1];
            const scalar z = pos[2];
            *pos++ = M(0,0)*x + M(0,1)*y + M(0,2)*z;
            *pos++ = M(1,0)*x + M(1,1)*y + M(1,2)*z;
            *pos++ = M(2,0)*x + M(2,1)*y + M(2,2)*z;
        }
#undef M
    }

/* from Desmond */
#define DET_3x3( m11, m12, m13,                                    \
        m21, m22, m23,                                    \
        m31, m32, m33 )                                   \
( ( (m11)*(m22)*(m33) + (m12)*(m23)*(m31) + (m13)*(m21)*(m32) ) \
  - ( (m11)*(m23)*(m32) + (m12)*(m21)*(m33) + (m13)*(m22)*(m31) ) )

    /* compute the alignment matrix by the method of Kabsch (1976, 1978).
     * Returned matrix is in row-major format, such that the matrix-vector
     * product R x_pos = x_ref; i.e. mat superposes pos onto ref.  ref 
     * and pos are assumed to be already centered on the origin.  There should
     * be n ref positions.  If wts is non-NULL, it should hold n non-negative
     * weights which sum to a nonzero value.  The alignment will be computed 
     * between ref * position i and pos position ids[i] */
    template <typename scalar, typename pos_scalar>
    double compute_alignment( 
            unsigned        n,      /* number of positions     */
            const unsigned *ids,    /* n ids                   */
            const scalar   *ref,    /* nx3 reference positions */
            const pos_scalar *pos,    /* max|ref| positions      */
            pos_scalar       *mat,    /* returned matrix         */
            const scalar   *wts=NULL  /* n optional weights */
            ) {


        double U[9], S[3], V[9], V_tr[9];
        double W[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
        double E0 = 0;
        unsigned i;
    
        for (i=0; i<n; i++) {
            const pos_scalar* p = ids ? pos+3*ids[i] : pos+3*i;
            const scalar w = wts ? *wts++ : 1;
            const scalar rx = ref[0];
            const scalar ry = ref[1];
            const scalar rz = ref[2];
            const pos_scalar px = p[0];
            const pos_scalar py = p[1];
            const pos_scalar pz = p[2];
            W[3*0 + 0] += w * rx * px;
            W[3*0 + 1] += w * rx * py;
            W[3*0 + 2] += w * rx * pz;
            W[3*1 + 0] += w * ry * px;
            W[3*1 + 1] += w * ry * py;
            W[3*1 + 2] += w * ry * pz;
            W[3*2 + 0] += w * rz * px;
            W[3*2 + 1] += w * rz * py;
            W[3*2 + 2] += w * rz * pz;
            E0 += w*rx*rx;
            E0 += w*ry*ry;
            E0 += w*rz*rz;
            E0 += w*px*px;
            E0 += w*py*py;
            E0 += w*pz*pz;
            ref += 3;
        }
    
        svd_3x3(W, S, V); 
        trans_3x3(V_tr, V);
    
        const double det_W = DET_3x3(W[0], W[1], W[2],
                W[3], W[4], W[5],
                W[6], W[7], W[8]);
    
        const double det_V_tr = DET_3x3(V_tr[0], V_tr[1], V_tr[2],
                V_tr[3], V_tr[4], V_tr[5],
                V_tr[6], V_tr[7], V_tr[8]);
    
        if(det_W * det_V_tr < 0) {
            W[3*0+2] *= -1;
            W[3*1+2] *= -1;
            W[3*2+2] *= -1;
            S[2]     *= -1;
        }
        matmult_3x3(U, W, V_tr);
    
        for (i=0; i<9; i++) mat[i] = (pos_scalar)U[i];
        return std::sqrt(std::fabs(E0-2*(S[0]+S[1]+S[2]))/n);
    }

    /* compute rmsd between ref and pos.  Arguments have the same semantics as
     * in pfx_compute_alignment.  center is added to each ref position */
    template <typename scalar, typename pos_scalar>
    scalar compute_rmsd(unsigned n,
                        const unsigned *ids,
                        const scalar *ref,
                        const scalar *center,
                        const pos_scalar *pos,
                        const scalar *wts=NULL) {
        if (n==0) return 0;
        scalar R=0, W=0;
        const scalar cx = center[0];
        const scalar cy = center[1];
        const scalar cz = center[2];
        for (; n; --n, ++ids, ref += 3) {
            const scalar w = wts ? *wts++ : 1;
            const pos_scalar *p = pos + 3*ids[0];
            scalar dx = (cx + ref[0]) - p[0];
            scalar dy = (cy + ref[1]) - p[1];
            scalar dz = (cz + ref[2]) - p[2];
            scalar r2 = dx*dx + dy*dy + dz*dz;
            R += w * r2;
            W += w;
        }
        return W ? std::sqrt(R/W) : 0;
    }

}}}

#endif

