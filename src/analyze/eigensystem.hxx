#ifndef desres_msys_analyze_eigensystem_hxx
#define desres_msys_analyze_eigensystem_hxx

#define COPY(dest,src,n) do {                                           \
    if ((dest)!=(src)) memcpy( (dest), (src), sizeof(*(src))*(n) );     \
  } while(0)


namespace desres { namespace msys {

    typedef struct tolerance {
      double absolute;
      double relative;
      int itmin;
      int itmax;
    } tolerance_t;
    
    void
    real_symmetric_eigenvalues_3x3( const double      * a0,
    				double            * opt_d,
    				double            * opt_v,
    				const tolerance_t * opt_tol );
}}

#endif 
