#include "../amber.hxx"
#include "molfile_homebox.h"
#include <fstream>
#include <boost/lexical_cast.hpp>

#include "molfile_plugin.h"

using namespace desres::msys;

static int open_rst(std::istream& in) {

    std::string line;

    /* read the title */
    if (!std::getline(in, line)) 
        MSYS_FAIL("Reading title of rst file");
    
    /* read the number of atoms.  Use atoi so that we can read rst7 too */
    if (!std::getline(in, line))
        MSYS_FAIL("Reading number of atoms of rst file");
    return atoi(line.c_str());
}

template <typename Float>
static void read_frame(std::istream& in, int natoms, 
                       Float* pos, Float* vel, double* box) {

    std::string line;
    /* looks to me like we've allotted 12 columns per float */
    static const int width = 12;
    for (int i=0; i<natoms;) {
        if (!std::getline(in, line))
            MSYS_FAIL("Reading coordinates of rst file");
        for (int j=0; j<2; j++) { /* at most two positions per line */
            try { 
                for (int k=0; k<3; k++) {
                    std::string s = line.substr(width*(3*j+k),width);
                    trim(s);
                    if (pos) *pos++ = boost::lexical_cast<Float>(s);
                }
                if (++i == natoms) break;
            }
            catch (std::exception& e) {
                MSYS_FAIL("Parsing positions in rst file: " << line
                        << "\n" << e.what());
            }
        }
    }
    /* velocities */
    static const Float SCALE = 20.455;
    for (int i=0; i<natoms;) {
        if (!std::getline(in, line)) {
            if (i==0) return;
            MSYS_FAIL("Reading velocities of rst file");
        }
        for (int j=0; j<2; j++) { /* at most two positions per line */
            try { 
                for (int k=0; k<3; k++) {
                    std::string s = line.substr(width*(3*j+k),width);
                    trim(s);
                    if (vel) *vel++ = SCALE * boost::lexical_cast<Float>(s);
                }
                if (++i == natoms) break;
            }
            catch (std::exception& e) {
                MSYS_FAIL("Parsing velocities in rst file: " << line
                        << "\n" << e.what());
            }
        }
    }
    /* unit cell */
    if (std::getline(in, line)) {
        /* read a, b, c, alpha, beta, gamma */
        Float pdbbox[6];
        for (int k=0; k<6; k++) {
            std::string s = line.substr(width*k,width);
            trim(s);
            try {
                pdbbox[k] = boost::lexical_cast<Float>(s);
            }
            catch (std::exception& e) {
                MSYS_FAIL("Parsing unit cell in rst file: " << line << "\n" << e.what());
            }
        }
        if (box) molfile_unitcell_from_pdb(box,
                pdbbox[0], pdbbox[1], pdbbox[2],
                pdbbox[3], pdbbox[4], pdbbox[5]);
    }
}

void desres::msys::ImportCrdCoordinates( SystemPtr mol, 
                                         std::string const& path ) {
    std::ifstream in(path.c_str());
    if (!in) MSYS_FAIL("Could not open rst file at " << path);

    Id natoms = open_rst(in);

    /* number must match molecule */
    if (natoms != mol->atomCount())
        MSYS_FAIL("Number of atoms " << natoms << " in rst file " << path << " does not match number in system " << mol->atomCount());

    std::vector<double> pos(3*natoms), vel(3*natoms);
    read_frame(in, natoms, &pos[0], &vel[0], mol->global_cell[0]);
    mol->setPositions(pos.begin());
    mol->setVelocities(vel.begin());
}    

static void* open_read(const char* path, const char* filetype,
                       int* natoms) {
    std::ifstream* in = new std::ifstream(path);
    if (!(*in)) {
        delete in;
        return NULL;
    }
    try {
        *natoms = open_rst(*in);
    }
    catch (std::exception& e) {
        fprintf(stderr, "%s\n", e.what());
        delete in;
        return NULL;
    }
    return in;
}

static void close_read(void *v) {
    delete (std::ifstream*)v;
}

static int read_next_timestep(void* v, int natoms, molfile_timestep_t* ts) {
    std::ifstream* in = (std::ifstream*)v;
    try {
        read_frame(*in, natoms, ts->coords, ts->velocities, ts->unit_cell);
    } 
    catch (std::exception& e) {
        fprintf(stderr, "%s\n", e.what());
        return MOLFILE_ERROR;
    }
    return MOLFILE_SUCCESS;
}

static int 
read_timestep_metadata(void *v, molfile_timestep_metadata *m) {
  m->has_velocities = 1;
  m->count = -1;
  m->supports_double_precision = false;
  return MOLFILE_SUCCESS;
}

static molfile_plugin_t plugin;
extern "C" int msys_rstplugin_init() {
    plugin.abiversion = vmdplugin_ABIVERSION;
    plugin.type = MOLFILE_PLUGIN_TYPE;
    plugin.name = "rst";
    plugin.prettyname = "Amber restart file";
    plugin.filename_extension = "rst,rst7";
    plugin.author = "D.E. Shaw Research";
    plugin.majorv = 1;
    plugin.minorv = 0;
    plugin.is_reentrant = VMDPLUGIN_THREADSAFE;
    plugin.open_file_read = open_read;
    plugin.read_timestep_metadata = read_timestep_metadata;
    plugin.read_next_timestep = read_next_timestep;
    plugin.close_file_read = close_read;
    return VMDPLUGIN_SUCCESS;
}

extern "C" int msys_rstplugin_register(void* v, vmdplugin_register_cb cb) {
    (*cb)(v, (vmdplugin_t *)&plugin);
    return VMDPLUGIN_SUCCESS;
}

