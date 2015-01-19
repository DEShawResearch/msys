#include "molfile_plugin.h"
#include "libmolfile_plugin.h"
#include "../io.hxx"

#include <boost/algorithm/string.hpp>
#include <stdio.h>

using namespace desres::msys;

/*
 * Implementation of molfile_plugin_t using msys::Load() 
 */

struct system_t {
    SystemPtr mol;
    std::vector<int> bonds;
    std::vector<float> bondorders;
    ssize_t timestep;
    explicit system_t(SystemPtr m) : mol(m), timestep() {}
};

static void* open_file_read(const char *filename, const char *filetype, 
        int *natoms ) {
    std::string uppercase_type(filetype);
    boost::to_upper(uppercase_type);
    FileFormat format = FileFormatFromString(uppercase_type);
    if (format == UnrecognizedFileFormat) {
        fprintf(stderr, "Unrecognized format '%s'\n", filetype);
        return NULL;
    }
    SystemPtr mol;
    const bool structure_only = true;
    try {
        mol = LoadWithFormat(filename, format, structure_only);
    }
    catch (std::exception& e) {
        fprintf(stderr, e.what());
        return NULL;
    }
    *natoms = mol->atomCount();
    return new system_t(mol);
}

static void close_file_read(void* v) {
    delete (system_t*)v;
}

static int read_structure(void *v, int *optflags, molfile_atom_t *atoms) {
    SystemPtr mol = ((system_t*)v)->mol;
    *optflags = MOLFILE_INSERTION 
              | MOLFILE_MASS
              | MOLFILE_CHARGE
              | MOLFILE_ATOMICNUMBER
              | MOLFILE_CTNUMBER
              ;

    for (Id i=0, n=mol->atomCount(); i<n; i++) {
        molfile_atom_t& dst = atoms[i];
        atom_t const& atm = mol->atomFAST(i);
        residue_t const& res = mol->residueFAST(atm.residue);
        chain_t const& chn = mol->chainFAST(res.chain);

        strncpy(dst.name, atm.name.c_str(), sizeof(dst.name));
        strcpy(dst.type, dst.name);
        strncpy(dst.resname, res.name.c_str(), sizeof(dst.resname));
        dst.resid = res.resid;
        strncpy(dst.segid, chn.segid.c_str(), sizeof(dst.segid));
        strncpy(dst.chain, chn.name.c_str(), sizeof(dst.chain));
        
        strncpy(dst.insertion, res.insertion.c_str(), sizeof(dst.insertion));
        dst.mass = atm.mass;
        dst.charge = atm.charge;
        dst.atomicnumber = atm.atomic_number;
        dst.ctnumber = chn.ct;

#define NULL_TERMINATE(x) do { x[sizeof(x)-1] = 0; } while(0)
        NULL_TERMINATE(dst.name);
        NULL_TERMINATE(dst.type);
        NULL_TERMINATE(dst.resname);
        NULL_TERMINATE(dst.segid);
        NULL_TERMINATE(dst.chain);
        NULL_TERMINATE(dst.altloc);
        NULL_TERMINATE(dst.insertion);
#undef  NULL_TERMINATE
    }
    return MOLFILE_SUCCESS;
}

static int read_bonds(void *v, int *nbonds, int **from, int **to, 
                      float **bondorder,
                      int **bondtype, int *nbondtypes, char ***bondtypename) {
    system_t* sys = (system_t*)v;
    SystemPtr mol = sys->mol;
    if (sys->bonds.empty()) {
        sys->bonds.resize(2*mol->bondCount());
        sys->bondorders.resize(mol->bondCount());
        for (Id i=0, n=mol->bondCount(); i<n; i++) {
            bond_t const& bnd = mol->bondFAST(i);
            sys->bonds[i]   = bnd.i+1;
            sys->bonds[i+n] = bnd.j+1;
            sys->bondorders[i] = bnd.resonant_order;
        }
    }
    *nbonds = mol->bondCount();
    if (*nbonds) {
        *from = &sys->bonds[0];
        *to   = &sys->bonds[*nbonds];
        *bondorder = &sys->bondorders[0];
    }
    return MOLFILE_SUCCESS;
}

static int read_next_timestep(void *v, int natoms, molfile_timestep_t *ts) {
    system_t* sys = (system_t*)v;
    if (sys->timestep++) return MOLFILE_EOF;
    SystemPtr mol = sys->mol;
    float* fpos = ts->coords;
    float* fvel = ts->velocities;
    double* dpos = ts->dcoords;
    double* dvel = ts->dvelocities;
    for (Id i=0, n=mol->atomCount(); i<n; i++) {
        Float* pos = &mol->atomFAST(i).x;
        Float* vel = &mol->atomFAST(i).vx;
        if (fpos) std::copy(pos, pos+3, fpos+3*i);
        if (fvel) std::copy(vel, vel+3, fvel+3*i);
        if (dpos) std::copy(pos, pos+3, dpos+3*i);
        if (dvel) std::copy(vel, vel+3, dvel+3*i);
    }
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

