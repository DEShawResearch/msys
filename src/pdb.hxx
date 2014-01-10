#ifndef desres_msys_pdb_hxx
#define desres_msys_pdb_hxx

#include "system.hxx"

namespace desres { namespace msys {

    SystemPtr ImportPDB( std::string const& path );
    SystemPtr ImportWebPDB( std::string const& code);

    void ImportPDBCoordinates( SystemPtr mol, std::string const& path );

    void ExportPDB(SystemPtr mol, std::string const& path);

    void ImportPDBUnitCell(double A, double B, double C,
                           double alpha, double beta, double gamma,
                           double* cell_3x3);

    void ExportPDBUnitCell(const double* cell_3x3,
                           double *A, double *B, double *C,
                           double *alpha, double *beta, double *gamma);

}}

#endif
