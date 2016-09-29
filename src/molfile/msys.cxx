#include "molfile_plugin.h"
#include "libmolfile_plugin.h"
#include "../io.hxx"
#include "../import.hxx"

#include <stdio.h>

using namespace desres::msys;

/*
 * Implementation of molfile_plugin_t using msys::Load() 
 */

struct system_t {
    SystemPtr mol;
    int natoms;
    std::string path;
    FileFormat format;

    std::vector<int> bonds;
    std::vector<float> bondorders;
    ssize_t timestep;
    explicit system_t(SystemPtr m) 
    : mol(m), natoms(), format(UnrecognizedFileFormat), timestep() {}
};

static void* open_file_read(const char *filename, const char *filetype, 
        int *natoms ) {
    std::string uppercase_type(filetype);
    to_upper(uppercase_type);
    FileFormat format = FileFormatFromString(uppercase_type);
    if (format == UnrecognizedFileFormat) {
        fprintf(stderr, "Unrecognized format '%s'\n", filetype);
        return NULL;
    }
    SystemPtr mol;
    const bool structure_only = false;
    const bool without_tables = true;
    try {
        mol = LoadWithFormat(filename, format, structure_only, without_tables);
    }
    catch (std::exception& e) {
        fprintf(stderr, "%s", e.what());
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

    Id occupancy = mol->atomPropIndex("occupancy");
    Id bfactor = mol->atomPropIndex("bfactor");
    Id altloc = mol->atomPropIndex("altloc");
    if (!bad(occupancy)) *optflags |= MOLFILE_OCCUPANCY;
    if (!bad(bfactor))   *optflags |= MOLFILE_BFACTOR;
    if (!bad(altloc))    *optflags |= MOLFILE_ALTLOC;

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

        if (!bad(occupancy)) dst.occupancy = mol->atomPropValue(i,occupancy);
        if (!bad(bfactor)) dst.bfactor = mol->atomPropValue(i,bfactor);
        if (!bad(altloc)) strncpy(dst.altloc, mol->atomPropValue(i,altloc).c_str(),sizeof(dst.altloc));

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
            sys->bondorders[i] = bnd.order;
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
    if (!ts) return MOLFILE_SUCCESS;
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
    memcpy(ts->unit_cell, mol->global_cell[0], sizeof(ts->unit_cell));
    return MOLFILE_SUCCESS;
}

static void* open_file_write(const char *filepath, const char *filetype,
          int natoms) {

    std::string uppercase_type(filetype);
    to_upper(uppercase_type);
    FileFormat format = FileFormatFromString(uppercase_type);
    if (format == UnrecognizedFileFormat) {
        fprintf(stderr, "Unrecognized format '%s'\n", filetype);
        return NULL;
    }
    SystemPtr mol = System::create();
    system_t* sys = new system_t(mol);
    sys->natoms = natoms;
    sys->format = format;
    sys->path = filepath;
    return sys;
}

static int write_structure(void *v, int optflags, const molfile_atom_t *atoms) {
    system_t* sys = (system_t*)v;
    SystemPtr mol = sys->mol;
    SystemImporter imp(mol);
    Id bfactor = optflags & MOLFILE_BFACTOR ?
        mol->addAtomProp("bfactor", FloatType) : BadId;
    Id occupancy = optflags & MOLFILE_OCCUPANCY ?
        mol->addAtomProp("occupancy", FloatType) : BadId;
    Id altloc = optflags & MOLFILE_ALTLOC?
        mol->addAtomProp("altloc", StringType) : BadId;

    for (int i=0, n=sys->natoms; i<n; i++) {
        molfile_atom_t const& atom = atoms[i];
        Id id = imp.addAtom(atom.chain, atom.segid,
                            atom.resid, atom.resname, atom.name,
                            atom.insertion, atom.ctnumber);
        atom_t& atm = mol->atomFAST(id);
        residue_t& res = mol->residueFAST(atm.residue);
        atm.mass = atom.mass;
        atm.charge = atom.charge;
        atm.atomic_number = atom.atomicnumber;
        res.insertion = atom.insertion;
        if (!bad(bfactor))   mol->atomPropValue(i,bfactor)=atom.bfactor;
        if (!bad(occupancy)) mol->atomPropValue(i,occupancy)=atom.occupancy;
        if (!bad(altloc))    mol->atomPropValue(i,altloc)=atom.altloc;
    }
    for (Id i=0, n=sys->bondorders.size(); i<n; i++) {
        Id id = mol->addBond(sys->bonds[i], sys->bonds[i+n]);
        mol->bondFAST(id).order = sys->bondorders[i];
    }
    return MOLFILE_SUCCESS;
}

static int write_timestep(void *v, const molfile_timestep_t *ts) {
    system_t* sys = (system_t*)v;
    SystemPtr mol = sys->mol;
    const float* fpos = ts->coords;
    const float* fvel = ts->velocities;
    const double* dpos = ts->dcoords;
    const double* dvel = ts->dvelocities;
    for (int i=0, n=sys->natoms; i<n; i++) {
        atom_t& atm = mol->atomFAST(i);
        Float* pos = &atm.x;
        Float* vel = &atm.vx;
        if (fpos)       std::copy(fpos+3*i, fpos+3*(i+1), pos);
        else if (dpos)  std::copy(dpos+3*i, dpos+3*(i+1), pos);
        if (fvel)       std::copy(fvel+3*i, fvel+3*(i+1), vel);
        else if (dvel)  std::copy(dvel+3*i, dvel+3*(i+1), vel);
    }
    try {
        static char* argv[] = {(char *)"vmd", 0};
        unsigned options = SaveOptions::Default;
        if (sys->timestep++) {
            options |= SaveOptions::Append;
        }
        SaveWithFormat(mol, sys->path, Provenance::fromArgs(1,argv),
                            sys->format, options);
    }
    catch (std::exception& e) {
        fprintf(stderr, "%s", e.what());
        return MOLFILE_ERROR;
    }
    return MOLFILE_SUCCESS;
}

static void close_file_write(void *v) {
    delete (system_t *)v;
}

static int write_bonds(void *v, 
                       int nbonds, int *from, int *to, float *bondorder,
                       int *bondtype, int nbondtypes, char **bondtypename) {
    system_t* sys = (system_t*)v;
    sys->bondorders.resize(nbonds);
    sys->bonds.resize(2*nbonds);
    for (int i=0, n=nbonds; i<nbonds; i++) {
        sys->bonds[i]   = from[i]-1;
        sys->bonds[n+i] =   to[i]-1;
        sys->bondorders[i] = bondorder[i];
    }
    return MOLFILE_SUCCESS;
}

static int read_timestep_metadata(void *v, molfile_timestep_metadata_t *meta) {
    system_t* sys = (system_t*)v;
    SystemPtr mol = sys->mol;
    meta->count = 1;
    meta->avg_bytes_per_timestep = mol->atomCount() * 48; /* pos+vel in dbl */
    meta->has_velocities = 1;
    meta->supports_double_precision = 1;
    return MOLFILE_SUCCESS;
}

static int read_timestep2(void* v, molfile_ssize_t index, 
                                   molfile_timestep_t* ts) {
    system_t* sys = (system_t*)v;
    SystemPtr mol = sys->mol;
    if (index!=0) return MOLFILE_EOF;
    ssize_t step = sys->timestep;
    sys->timestep = 0;
    int rc = read_next_timestep(v, mol->atomCount(), ts);
    sys->timestep = step;
    return rc;
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
     close_file_read,
     open_file_write,
     write_structure,
     write_timestep,
     close_file_write,
     NULL, // sync_file_write,
     NULL, // read_volumetric_metadata,
     NULL, // read_volumetric_data,
     NULL, // read_rawgraphics,
     NULL, // read_molecule_metadata,
     write_bonds,
     NULL, // write_volumetric_data
     NULL, // read_angles
     NULL, // write_angles
     NULL, // read_qm_metadata
     NULL, // read_qm_rundata
     NULL, // read_timestep
     read_timestep_metadata,
     NULL, // read_qm_timestep_metadata,
     read_timestep2,
     NULL, // read_times
     NULL, // cons_fputs
     NULL  // truncate_file_write
};

static void init(molfile_plugin_t* self, 
                 const char* name, const char* pname, const char* ext) {
    *self = msys_plugin_base;
    self->name = name;
    self->prettyname = pname;
    self->filename_extension = ext;
}

static molfile_plugin_t dmsplugin[1];
static molfile_plugin_t maeplugin[1];
static molfile_plugin_t pdbplugin[1];
static molfile_plugin_t webpdbplugin[1];
static molfile_plugin_t prm7plugin[1];
static molfile_plugin_t mol2plugin[1];
static molfile_plugin_t sdfplugin[1];

extern "C"
int msys_plugin_init() {
    init(dmsplugin, "dms",  "DESRES Structure", "dms,dms.gz");
    init(maeplugin, "mae",  "MAESTRO file",     "mae,mae.gz,maegz,maeff,maeff.gz,cms,cms.gz");
    init(pdbplugin, "pdb",  "PDB",              "pdb");
    init(webpdbplugin,"webpdb","Web PDB Download", "");
    init(prm7plugin,"parm7","AMBER7 Parm",      "prmtop,prm7,parm7");
    init(mol2plugin,"mol2", "MDL mol2",         "mol2");
    init(sdfplugin, "sdf",  "SDF",              "sdf,sdf.gz,sdfgz");
    return VMDPLUGIN_SUCCESS;
}

extern "C"
int msys_plugin_register(void* v, vmdplugin_register_cb cb) {
      cb( v, (vmdplugin_t *)dmsplugin);
      cb( v, (vmdplugin_t *)maeplugin);
      cb( v, (vmdplugin_t *)pdbplugin);
      cb( v, (vmdplugin_t *)webpdbplugin);
      cb( v, (vmdplugin_t *)prm7plugin);
      cb( v, (vmdplugin_t *)mol2plugin);
      cb( v, (vmdplugin_t *)sdfplugin);
      return VMDPLUGIN_SUCCESS;
}

