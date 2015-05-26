#ifndef desres_msys_amber_hxx
#define desres_msys_amber_hxx

#include "system.hxx"

namespace desres { 
namespace msys {

    SystemPtr ImportPrmTop( std::string const& path,
                            bool structure_only,
                            bool without_tables);

    /* Create a new System from an Amber7 prmtop file.  The system's 
       coordinates and box will be all zero. */
    inline
    SystemPtr ImportPrmTop( std::string const& path,
                            bool structure_only = false) {
        return ImportPrmTop(path, structure_only, structure_only);
    }
    
    /* Read coordinates from a crd file into an existing System, which
       must have the same number of atoms as the file. */
    void ImportCrdCoordinates( SystemPtr mol, std::string const& path );

}}

#endif
