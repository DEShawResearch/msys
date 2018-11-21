/* @COPYRIGHT@ */

#include "geom.hxx"
#include "analyze/eigensystem.hxx"
#include <math.h>
#include <stdio.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795029
#endif

using namespace desres::msys;

static double calc_dot_prod( const double* A, const double* B ) {
    return A[0]*B[0] + A[1]*B[1] + A[2]*B[2];
}

namespace desres { namespace msys {

    double calc_distance( const double* A, const double* B ) {
        double dx=A[0]-B[0];
        double dy=A[1]-B[1];
        double dz=A[2]-B[2];
        return sqrt(dx*dx + dy*dy + dz*dz);
    }
    
    void calc_cross_prod( double* C, const double* A, const double* B ) {
        double cx, cy, cz;
        cx =  A[1]*B[2] - B[1]*A[2];
        cy = -A[0]*B[2] + B[0]*A[2];
        cz =  A[0]*B[1] - B[0]*A[1];
        C[0] = cx;
        C[1] = cy;
        C[2] = cz;
    }
    
    double calc_normal( double* A ) {
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
    
    double calc_vec_angle( const double* r1, const double* r2 ) {
        double psin, pcos;
        double r3[3];
        calc_cross_prod(r3, r1, r2 );
        psin = sqrt(calc_dot_prod(r3,r3));
        pcos = calc_dot_prod( r1, r2 );
        return atan2(psin,pcos);
    }
    
    double calc_angle( const double* A, const double* B, const double* C ) {
        int i;
        double r1[3], r2[3];
        for (i=0; i<3; i++) {
            r1[i]=A[i]-B[i];
            r2[i]=C[i]-B[i];
        }
        return calc_vec_angle( r1, r2 );
    }

    double calc_vec_dihedral( const double* r1, const double* r2,
                              const double* r3 ) {
        double n1[3], n2[3];
        double psin, pcos;
        calc_cross_prod(n1, r1, r2 );
        calc_cross_prod(n2, r2, r3 );
        psin = calc_dot_prod( n1, r3 ) * sqrt(calc_dot_prod( r2, r2 ));
        pcos = calc_dot_prod( n1, n2 );
        return atan2(psin,pcos);
    }
        
    double calc_dihedral( const double* A, const double* B, 
                          const double* C, const double* D ) {
        int i;
        double r1[3],r2[3],r3[3];
        for (i=0; i<3; i++) {
            r1[i]=B[i]-A[i];
            r2[i]=C[i]-B[i];
            r3[i]=D[i]-C[i];
        }
        return calc_vec_dihedral(r1, r2, r3);
    }
    
    void apply_dihedral_geometry( double* D, const double* tmpA,
                                  const double* B,
                                  const double* C,
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
        double a[3],b[3], t[3], u[3];
        double A[3] = {tmpA[0], tmpA[1], tmpA[2]};
        double bn=0, tn=0, un=0;
        double alpha, beta, gamma;
        int i;
    
        if (B==C) {
            D[0]=C[0]+r;
            D[1]=C[1];
            D[2]=C[2];
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
        calc_cross_prod( t, a, b );
        calc_cross_prod( u, b, t );
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
    }

    double calc_planarity(int n, const double* pos) {
        if (n<3) return 0;
        const double* end = pos+3*n;
 
        /* compute center */
        double xc=0, yc=0, zc=0;
        for (const double* p = pos; p!=end; p+=3) {
            xc += p[0];
            yc += p[1];
            zc += p[2];
        }
        xc /= n;
        yc /= n;
        zc /= n;

        /* compute inertia tensor */
        double I[9]={0,0,0,0,0,0,0,0,0};
        for (const double* p = pos; p!=end; p+=3) {
            double x = p[0] - xc;
            double y = p[1] - yc;
            double z = p[2] - zc;
            I[0] += (y*y + z*z);
            I[1] += -(x*y);
            I[2] += -(x*z);
            I[4] += (x*x + z*z);
            I[5] += -(y*z);
            I[8] += (x*x + y*y);
        }
        I[3]=I[1];
        I[6]=I[2];
        I[7]=I[5];
        
        double v[3];
        real_symmetric_eigenvalues_3x3(I,v,NULL,NULL);
        // perpendicular axis theorem
        double planar=fabs(v[0]-(v[1]+v[2]));
        return planar;
    }

    bool line_intersects_tri( const double* A,
                              const double* B,
                              const double* C,
                              const double* R,
                              const double* S) {

        double RS[3], AB[3], BC[3], norm[3], LT[3], I[3];
        for (int i=0; i<3; i++) {
            RS[i] = S[i]-R[i];
            AB[i] = B[i]-A[i];
            BC[i] = C[i]-B[i];
        }
        calc_cross_prod(norm, AB, BC);
        double norm_dot_line = calc_dot_prod(norm, RS);
        if (norm_dot_line==0) {
            /* line is parallel to triangle */
            return false;
        }

        for (int i=0; i<3; i++) LT[i] = R[i]-A[i];
        double t = -calc_dot_prod(norm, LT) / norm_dot_line;
        if (t<=0 || t>=1) {
            /* intersection is past line's endpoints */
            return false;
        }

        for (int i=0; i<3; i++) I[i] = R[i] + RS[i]*t;

        double areas[3];
        for (int i=0; i<3; i++) {
            const double* t0 = i==0 ? A : i==1 ? B : C;
            const double* t1 = i==0 ? B : i==1 ? C : A;
            double v0[3], v1[3], c[3];
            for (int j=0; j<3; j++) {
                v0[j] = t1[j]-t0[j];
                v1[j] = I[j] - t0[j];
            }
            calc_cross_prod(c, v0, v1);
            areas[i] = calc_dot_prod(c, norm);

            if (areas[i]==0) return false;
        }
        return signbit(areas[0])==signbit(areas[1]) &&
               signbit(areas[0])==signbit(areas[2]);
    }

    void calc_projection(const double* c1, const double* c2,
                         const double* c3, const double* p,
                         double* proj) {

        /* compute a normal */
        double r1[3], r3[3], p1[3];
        for (int i=0; i<3; i++) {
            r1[i] = c1[i] - c2[i];
            r3[i] = c3[i] - c2[i];
            p1[i] = c1[i] -  p[i];
        }

        calc_cross_prod(r1,r3,proj);
        calc_normal(proj);

        double s = calc_dot_prod(proj, p1);
        for (int i=0; i<3; i++) {
            proj[i] = proj[i] * s + p[i];
        }
    }


}}

