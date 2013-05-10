/* @COPYRIGHT@ */

namespace desres { namespace msys {

    double calc_distance( const double* A, const double* B );
    double calc_normal( double* A );
    
    /* Calculate angle between vectors in radians */
    double calc_vec_angle( const double* r1, const double* r2 );
    
    /* Calculate angle between positions in radians */
    double calc_angle( const double* A, const double* B, const double* C );
    
    /* Calculate dihedral between postions in radians */
    double calc_dihedral( const double* A, const double* B, 
                          const double* C, const double* D );
    
    /* Compute cross product; store in dst */
    void calc_cross_prod( double* dst, const double* A, const double* B);

    /* Reconstruct the position of D from positions A, B, and C, and
     * CD length r, BCD angle theta, ABCD dihedral phi.  */
    void apply_dihedral_geometry( double* dst,
                                  const double* A,
                                  const double* B,
                                  const double* C,
                                  double r, double theta, double phi );

}}

