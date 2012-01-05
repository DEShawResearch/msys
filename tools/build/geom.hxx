/* @COPYRIGHT@ */

#include <math.h>

#include "../../src/system.hxx"

#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795029
#endif

namespace desres { namespace msys {

    double calc_distance( const Vec3& A, const Vec3& B );
    Vec3 calc_cross_prod( const Vec3& A, const Vec3& B );
    double calc_normal( Vec3& A );
    
    /* Calculate angle between vectors in radians */
    double calc_vec_angle( const Vec3& r1, const Vec3& r2 );
    
    /* Calculate angle between positions in radians */
    double calc_angle( const Vec3& A, const Vec3& B, const Vec3& C );
    
    /* Calculate dihedral between postions in radians */
    double calc_dihedral( const Vec3& A, const Vec3& B, 
                          const Vec3& C, const Vec3& D );
    
    /* Reconstruct the position of D from positions A, B, and C, and
     * CD length r, BCD angle theta, ABCD dihedral phi.  */
    Vec3 apply_dihedral_geometry( const Vec3& A,
                                        const Vec3& B,
                                        const Vec3& C,
                                        double r, double theta, double phi );

}}

