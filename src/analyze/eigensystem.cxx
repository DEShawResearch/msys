#include <cmath>
#include <stdexcept>
#include <cstring>

#include "eigensystem.hxx" 

/* Given a real symmetric 3x3 matrix (stored uncompactly in either C
   or FORTRAN ordering), compute a diagonal matrix D and an
   orthonormal matrix V such that:

     A V = V D

   Thus, the k-th diagonal element of D is an eigenvalue of A and the
   k-th column of V is the corresponding eigenvector.

   If d is non-NULL, this routine will save the diagonal elements of
   the matrix D in d in descending order. (Note: If A is also positive
   definite, this corresponds to organizing by magnitude.)

   If v is non-NULL, this routine will save the matrix V in FORTRAN
   indexed order in v. (Thus, v[0:2] are the first eigenvector,
   v[3:5] are the second and v[6:8] are the third.)

   If tol is non-NULL, this routine will use a user specified
   convergence criterion. The convergence criterion is that the sum of
   the off-diagonal components of:

     || off( V^T A V ) ||_F <= relative_tolerance ||A||_F

   where off(M) means matrix formed by zeroing the diagonal of M and
   ||M||_F is the Frobenius norm of the matrix M (the sqrt of the sum
   of the squares of the matrix components). The default settings are
   to try to converge to 1e-7 in less then 100 Jacobi transformations.

   This routine very numerically robust to the point that it is
   _exactly_ permutation symmetric (convergence testing is not
   fortified against internal overflows). It is based on the classical
   Jacobi iteration. While classical Jacobi is highly inefficient for
   large matrices, it is acceptable here (the search process of
   classical Jacobi is what gives this routine its exact permutation
   symmetry). */
   
namespace desres { namespace msys {
void
real_symmetric_eigenvalues_3x3( const double      * a0,
                                double            * d,
                                double            * v,
                                const tolerance_t * tol ) {

# define a(i,j) a[(i)+3*(j)]
# define d(i)   d[(i)]
# define v(i,j) v[(i)+3*(j)]

  static const double identity_3x3[9] = { 1, 0, 0,
                                          0, 1, 0,
                                          0, 0, 1 };
  static const tolerance_t default_tol[1] = { { 0, 1e-7, 0, 100 } };

  double a[9];       /* A = V^T A0 V */
  double stack_d[3]; /* Storage in case d not user provided */
  double stack_v[9]; /* Storage in case v not user provided */
  double f0, f1;     /* Temporaries */
  double app, aqq, apq, aqr, arp; /* Components undergoing Jacobi */
  double atol = 0;    /* Absolute convergence criteria */
  double s, t, th;   /* sin theta, tan theta, tan theta/2 of Jacobi */
  int p, q, r;       /* p and q are Jacobi rows, r is other row */
  int k = 0;         /* Number of rotations performed so far */

  if( d==NULL && v==NULL ) return ; /* Nothing to do */

  if( a0==NULL ) throw std::runtime_error( "Bad input args");
  COPY( a, a0, 9 );
  if( a(0,1)!=a(1,0) ||
      a(1,2)!=a(2,1) ||
      a(2,0)!=a(0,2) )
      throw std::runtime_error("Matrix is not symmetric");

  if( d==NULL   ) d   = stack_d;
  if( v==NULL   ) v   = stack_v;
  if( tol==NULL ) tol = default_tol;       /* Use default tolerances */

  COPY( v, identity_3x3, 9 );

  for(;;) {

    /* Find the maximum magnitude off-diagonal element. */
    
    f1 = fabs( a(0,1) );             p = 0, q = 1, r = 2, f0 = f1;
    f1 = fabs( a(1,2) ); if( f1>f0 ) p = 1, q = 2, r = 0, f0 = f1;
    f1 = fabs( a(2,0) ); if( f1>f0 ) p = 2, q = 0, r = 1, f0 = f1;

    app = a(p,p);
    aqq = a(q,q);
    apq = a(p,q);
    aqr = a(q,r);
    arp = a(r,p);

    /* Check if we have converged or have performed the maximum number of
       Jacobi transformations allowed. Note that below:
         f0   = 0.5 || off(V^T A0 V) ||_F^2  and
         atol = 0.5 rtol^2 || A0 ||_F^2
       Thus f0<=atol is equivalent to the convergence criterion:
         || off(V^T A0 V) ||_F < rtol || A0 ||_F */

    f0 = apq*apq + aqr*aqr + arp*arp;
    if( k==0 )
      atol = tol->relative * tol->relative *
        ( 0.5*( app*app + aqq*aqq + a(r,r)*a(r,r) ) + f0 );
    if( f0<=atol || k==tol->itmax ) break;

    /* Compute the Jacobi transformation J = J(p,q,theta). Note that
       apq is guaranteed to be non-zero if we have not converged. */

    th = ( aqq - app ) / ( apq + apq ); /* th = cot 2 theta */
    s  = sqrt( 1 + th*th );
    if( th<0 ) s = -s;
    t  = 1 / ( th + s );                /* t  = tan theta */
    th = 1 / sqrt( 1 + t*t );           /* th = cos theta */
    s  = t*th;                          /* s  = sin theta */
    th = s / ( 1 + th );                /* th = tan theta/2 */
    
    /* Apply the Jacobi transformation to A ( A = J^T A J ) */

    a(p,p) = app - t*apq;
    a(q,q) = aqq + t*apq;
    a(p,q) = 0;                        /* a(q,p) = 0;      */
    a(q,r) = aqr + s*( arp - th*aqr ); /* a(r,q) = a(q,r); */
    a(r,p) = arp - s*( aqr + th*arp ); /* a(p,r) = a(r,p); */

    /* Accumulate the Jacobi transformation to V ( V = V J ) */

    f0 = v(0,p);
    f1 = v(0,q);
    v(0,p) = f0 - s*( f1 + th*f0 );
    v(0,q) = f1 + s*( f0 - th*f1 );

    f0 = v(1,p);
    f1 = v(1,q);
    v(1,p) = f0 - s*( f1 + th*f0 );
    v(1,q) = f1 + s*( f0 - th*f1 );

    f0 = v(2,p);
    f1 = v(2,q);
    v(2,p) = f0 - s*( f1 + th*f0 );
    v(2,q) = f1 + s*( f0 - th*f1 );

    /* Increment the Jacobi transformation counter */

    k++;
  }

  /* Extract the eigenvalues and put them in descending order */
  /* Note: At this point, f0<=atol can be used to test for convergence */

# define SWAP_ROWS(p,q)                                                 \
    if( p!=q ) {                                                        \
      /**/         d(q)   = d(p);   d(p)   = f1;                        \
      f1 = v(0,q); v(0,q) = v(0,p); v(0,p) = f1;                        \
      f1 = v(1,q); v(1,q) = v(1,p); v(1,p) = f1;                        \
      f1 = v(2,q); v(2,q) = v(2,p); v(2,p) = f1;                        \
    }

  d(0) = a(0,0);
  d(1) = a(1,1);
  d(2) = a(2,2);

  /**/          q = 0, f1 = d(0);
  if( d(1)>f1 ) q = 1, f1 = d(1);
  if( d(2)>f1 ) q = 2, f1 = d(2);
  SWAP_ROWS(0,q);

  /**/          q = 1, f1 = d(1);
  if( d(2)>f1 ) q = 2, f1 = d(2);
  SWAP_ROWS(1,q);

  if (f0>atol) throw std::runtime_error( "Did not converge to tolerance" );
  return ;

# undef a
# undef d
# undef v
# undef SWAP_ROWS

}
}} 



