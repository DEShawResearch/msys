/* @COPYRIGHT@ */

#include "geom.hxx"
#include <math.h>
#include <stdio.h>
#include <string.h>

using namespace desres::msys;

static double calc_dot_prod( const Vec3& A, const Vec3& B ) {
    return A[0]*B[0] + A[1]*B[1] + A[2]*B[2];
}

namespace desres { namespace msys {

    double calc_distance( const Vec3& A, const Vec3& B ) {
        double dx=A[0]-B[0];
        double dy=A[1]-B[1];
        double dz=A[2]-B[2];
        return sqrt(dx*dx + dy*dy + dz*dz);
    }
    
    Vec3 calc_cross_prod( const Vec3& A, const Vec3& B ) {
        Vec3 C;
        double cx, cy, cz;
        cx =  A[1]*B[2] - B[1]*A[2];
        cy = -A[0]*B[2] + B[0]*A[2];
        cz =  A[0]*B[1] - B[0]*A[1];
        C[0] = cx;
        C[1] = cy;
        C[2] = cz;
        return C;
    }
    
    double calc_normal( Vec3& A ) {
        double x=A[0];
        double y=A[1];
        double z=A[2];
        double norm=sqrt(x*x + y*y + z*z);
        if (norm>0) {
            A[0]=x/norm;
            A[1]=y/norm;
            A[2]=z/norm;
        }
        return norm;
    }
    
    double calc_vec_angle( const Vec3& r1, const Vec3& r2 ) {
        double psin, pcos;
    
        Vec3 r3 = calc_cross_prod( r1, r2 );
        psin = sqrt(calc_dot_prod(r3,r3));
        pcos = calc_dot_prod( r1, r2 );
        return atan2(psin,pcos);
    }
    
    double calc_angle( const Vec3& A, const Vec3& B, const Vec3& C ) {
        int i;
        Vec3 r1, r2;
        for (i=0; i<3; i++) {
            r1[i]=A[i]-B[i];
            r2[i]=C[i]-B[i];
        }
        return calc_vec_angle( r1, r2 );
    }
    
    double calc_dihedral( const Vec3& A, const Vec3& B, 
                          const Vec3& C, const Vec3& D ) {
        int i;
        Vec3 r1,r2,r3;
        double psin, pcos;
        for (i=0; i<3; i++) {
            r1[i]=B[i]-A[i];
            r2[i]=C[i]-B[i];
            r3[i]=D[i]-C[i];
        }
        Vec3 n1 = calc_cross_prod( r1, r2 );
        Vec3 n2 = calc_cross_prod( r2, r3 );
        psin = calc_dot_prod( n1, r3 ) * sqrt(calc_dot_prod( r2, r2 ));
        pcos = calc_dot_prod( n1, n2 );
        return atan2(psin,pcos);
    }
    
    Vec3 apply_dihedral_geometry( const Vec3& tmpA,
                                        const Vec3& B,
                                        const Vec3& C,
                                        double r, double theta, double phi ) {
    
        /* Let a=B-A, b=C-B, c=D-C.  Define an orthogonal coordinate system
         * [b, t = a x b, u = b x t].  We can express 
         *      c = beta b_hat + alpha t_hat + gamma u_hat
         *  for constants beta, alpha, gamma in terms of r, theta, and phi.
         *
         *  After some algebra: 
         *      
         *      beta  = - r cos theta
         *      alpha =   r sin theta sin phi
         *      gamma = - r sin theta cos phi
         */
        Vec3 a,b;
        Vec3 A(tmpA);
        Vec3 D;
        double bn=0, tn=0, un=0;
        double alpha, beta, gamma;
        int i;
    
        if (B==C) {
            D[0]=C[0]+r;
            D[1]=C[1];
            D[2]=C[2];
            return D;
        }
        if (A==C) {
            /* pick A not along B-C */
            A[0] = B[0]+1;
            A[1] = B[1];
            A[2] = B[2];
            if (B[1]==C[1] && B[2]==C[2]) {
                /* B and C are along x_hat; fine, we'll use y_hat */
                A[0] = B[0];
                A[1] = B[1]+1;
            }
        }
    
        for (i=0; i<3; i++) {
            a[i] = B[i] - A[i];
            b[i] = C[i] - B[i];
        }
        Vec3 t = calc_cross_prod( a, b );
        Vec3 u = calc_cross_prod( b, t );
        for (i=0; i<3; i++) {
            bn += b[i] * b[i];
            tn += t[i] * t[i];
            un += u[i] * u[i];
        }
        bn=sqrt(bn);
        tn=sqrt(tn);
        un=sqrt(un);
        if (bn>0) for (i=0; i<3; i++) b[i] /= bn;
        if (tn>0) for (i=0; i<3; i++) t[i] /= tn;
        if (un>0) for (i=0; i<3; i++) u[i] /= un;
    
        beta  = -r * cos(theta);
        alpha =  r * sin(theta) * sin(phi);
        gamma = -r * sin(theta) * cos(phi);
    
        for (i=0; i<3; i++) {
            D[i] = C[i] + beta*b[i] + alpha*t[i] + gamma*u[i];
        }
        return D;
    }

}}
