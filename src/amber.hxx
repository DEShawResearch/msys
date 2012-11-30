#ifndef desres_msys_amber_hxx
#define desres_msys_amber_hxx

#include "system.hxx"

namespace desres { 
namespace msys {

    /* Create a new System from an Amber7 prmtop file.  The system's 
       coordinates and box will be all zero. */
    SystemPtr ImportPrmTop( std::string const& path,
                            bool structure_only = false);
    
    /* Read coordinates from a crd file into an existing System, which
       must have the same number of atoms as the file. */
    void ImportCrdCoordinates( SystemPtr mol, std::string const& path );

}}

#endif
