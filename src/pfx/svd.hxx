#ifndef desres_pfx_svd_hxx
#define desres_pfx_svd_hxx

#include <math.h>
#include <string.h>

namespace desres { namespace msys { namespace pfx {

    /* Swap columns p and q in 3x3 matrix A */
    static inline
    void swap_columns(double *A, int p, int q) {
        double x,y,z;
        x = A[p+0];
        y = A[p+3];
        z = A[p+6];
        A[p+0] = A[q+0];
        A[p+3] = A[q+3];
        A[p+6] = A[q+6];
        A[q+0] = x;
        A[q+3] = y;
        A[q+6] = z;
    }
    
    /* Find the Jacobi matrix [[c,s],[-s,c]] such that
     * J' A J  is diagonal, where A = [[x,y],[y,z]].  Set J = {c,s}. */
    static  inline
    void make_jacobi(double x, double y, double z, double *J) {
        if (y==0) {
            J[0] = 1;
            J[1] = 0;
        }  else {
            double Q = (x-z)/(2*fabs(y));
            double w = sqrt(1+Q*Q);
            double t = Q>0 ? 1/(Q+w) : 1/(Q-w);
            double c = 1/sqrt(1+t*t);
            double s = -copysign(1,y) * t * c;
    
            J[0] = c;
            J[1] = s;
        }
    }
    
    static inline
    void jacobi_svd(double a11, double a12, double a21, double a22,
            double *left_c_s, double *right_c_s) { 
    
        double t = a11 + a22;
        double d = a21 - a12;
        double c, s;
        double r11, r12, r22;
    
        if (t==0) {
            c = 0;
            s = copysign(1,d);
        } else {
            double u = d/t;
            c = 1/sqrt(1+u*u);
            s = c*u;
        }
    
        /* apply Jacobi (c,s) on the left of A to get r */
        r11 = a11*c + a21*s;
        r12 = a12*c + a22*s;
        //r21 = a21*c - a11*s;
        r22 = a22*c - a12*s;
    
        make_jacobi(r11,r12,r22, right_c_s);
    
        /* concatenate */
        left_c_s[0] = c*right_c_s[0] + s*right_c_s[1];
        left_c_s[1] = s*right_c_s[0] - c*right_c_s[1];
    }
    
    static inline
    void apply_jacobi_left(double *A, int p, int q, double c, double s) {
        int i;
        for (i=0; i<3; i++) {
            double x = A[3*p+i];
            double y = A[3*q+i];
            A[3*p+i] = s*y + c*x;
            A[3*q+i] = c*y - s*x;
        }
    }
    
    /* apply Jacobi on the right to columns p, q */
    static inline
    void apply_jacobi_right(double *A, int p, int q, double c, double s) {
        int i;
        for (i=0; i<3; i++) {
            double x = A[3*i+p];
            double y = A[3*i+q];
            A[3*i+p] = c*x - s*y;
            A[3*i+q] = s*x + c*y;
        }
    }
    
    /* machine epsilon */
    #ifndef DBL_EPSILON
    #define DBL_EPSILON 2.2204460492503131e-16
    #endif
    
    /* DBL_TRUE_MIN defined in C11 */
    #ifndef DBL_TRUE_MIN
    #define DBL_TRUE_MIN 4.9406564584124654e-324
    #endif
    
    inline
    void svd_3x3(double *A, double *w, double *V) {
    
        int done = 0;
        static const int diag = 3;
        static const double prec = 2*DBL_EPSILON;
        static const double zero = 2*DBL_TRUE_MIN;
        int i, j, p, q, kmax;
        double U[9]={1,0,0, 0,1,0, 0,0,1};
        memcpy(V,U,sizeof(U));
    
        while (!done) {
            done = 1;
    
            for (p=1; p<diag; p++) {
                for (q=0; q<p; q++) {
                    /* is this 2x2 sub-matrix diagonal? */
                    const double a = fabs(A[4*p]);
                    const double b = fabs(A[4*q]);
                    double threshold = prec * (a > b ? a : b);
                    /* small denormalized values treated as zero */
                    threshold = threshold > zero ? threshold : zero;
                    const double pq = fabs(A[3*p+q]);
                    const double qp = fabs(A[3*q+p]);
                    double offdiag = pq > qp ? pq : qp;
                    if (offdiag > threshold) {
                        done = 0;
                        double j_left[2], j_right[2], c,s;
                        jacobi_svd(A[4*p], A[3*p+q], A[3*q+p], A[4*q], 
                                j_left, j_right);
    
                        /* apply on the left, affecting rows p and q */
                        c = j_left[0];
                        s = j_left[1];
                        apply_jacobi_left(A,p,q,c,s);
                        apply_jacobi_right(U,p,q,c,-s);
    
                        /* apply on the right, affecting cols p and q */
                        c = j_right[0];
                        s = j_right[1];
                        apply_jacobi_right(V,p,q,c,s);
                        apply_jacobi_right(A,p,q,c,s);
                    }
                }
            }
        }
    
        /* make singular values positive, and find biggest */
        kmax = 0;
        for (i=0; i<3; i++) {
            double s = A[4*i];
            if (s<0) {
                s=-s;
                for (j=0; j<3; j++) U[3*j+i] *= -1;
            }
            if (i==0 || s>w[kmax]) kmax = i;
            w[i] = s;
        }
    
        /* sort */
        if (kmax!=0) {
            { double tmp = w[0]; w[0] = w[kmax]; w[kmax] = tmp; }
            swap_columns(U,0,kmax);
            swap_columns(V,0,kmax);
        }
        if (w[1]<w[2]) {
            { double tmp = w[1]; w[1] = w[2]; w[2] = tmp; }
            swap_columns(U,1,2);
            swap_columns(V,1,2);
        }
    
        memcpy(A,U,sizeof(U));
    }

}}}

#endif

