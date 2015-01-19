#include "molfile_plugin.h"
#include "libmolfile_plugin.h"
#include "../io.hxx"

/*
 * Implementation of molfile_plugin_t using msys::Load() 
 */

static void* open_file_read(const char *filename, const char *filetype, 
        int *natoms ) {
    return NULL;
}

static void close_file_read(void* v) {
}

static int read_structure(void *v, int *optflags, molfile_atom_t *atoms) {
    return MOLFILE_SUCCESS;
}

static int read_bonds(void *v, int *nbonds, int **from, int **to, 
                      float **bondorder,
                      int **bondtype, int *nbondtypes, char ***bondtypename) {
    return MOLFILE_SUCCESS;
}

static int read_next_timestep(void *v, int natoms, molfile_timestep_t *ts) {
    return MOLFILE_SUCCESS;
}

static molfile_plugin_t msys_plugin_base = {
     vmdplugin_ABIVERSION,
     MOLFILE_PLUGIN_TYPE,
     "give me a name",
     "give me a pretty name",
     "D.E. Shaw Research, LLC",
     1,
     0,
     VMDPLUGIN_THREADSAFE,
     "filename extensions",
     open_file_read,
     read_structure,
     read_bonds,
     read_next_timestep,
     close_file_read
};

struct msys_plugin_t : molfile_plugin_t {
    msys_plugin_t(const char* name, const char* pname, const char* ext) {
        static_cast<molfile_plugin_t&>(*this) = msys_plugin_base;
        this->name = name;
        this->prettyname = pname;
        this->filename_extension = ext;
    };
};

static msys_plugin_t dmsplugin( "dms",  "DESRES Structure", "dms,dms.gz");
static msys_plugin_t psfplugin( "psf",  "PSF",              "psf");
static msys_plugin_t mol2plugin("mol2", "MDL mol2",         "mol2");
static msys_plugin_t sdfplugin( "sdf",  "SDF",              "sdf,sdf.gz,sdfgz");
static msys_plugin_t xyzplugin( "xyz",  "XYZ",              "xyz");

int msys_plugin_register(void* v, vmdplugin_register_cb cb) {
      cb( v, (vmdplugin_t *)&dmsplugin);
      cb( v, (vmdplugin_t *)&psfplugin);
      cb( v, (vmdplugin_t *)&mol2plugin);
      cb( v, (vmdplugin_t *)&sdfplugin);
      cb( v, (vmdplugin_t *)&xyzplugin);
      return VMDPLUGIN_SUCCESS;
}

