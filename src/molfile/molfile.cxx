#include "molfile.hxx"
#include "findframe.hxx"
#include "libmolfile_plugin.h"

#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <limits.h>

#include <map>
#include <stdexcept>
#include <string>

#ifdef WIN32
#ifdef _WIN64
 typedef __int64 ssize_t;
#else
 typedef int ssize_t;
#endif

#endif


namespace ff = desres::molfile::findframe;

namespace {
    // mapping from plugin type to plugin instance
    std::map<std::string,molfile_plugin_t*>& plugindict() {
        static std::map<std::string,molfile_plugin_t*> dict;
        return dict;
    }
    // mapping from filename extension to plugin instance
    typedef std::multimap<std::string,molfile_plugin_t*> ExtensionDict;
    ExtensionDict& extensiondict() {
        static ExtensionDict dict;
        return dict;
    }

    int register_cb(void *v, vmdplugin_t *p) {
        molfile_plugin_t *plugin = reinterpret_cast<molfile_plugin_t*>(p);
        plugindict()[plugin->name] = plugin;
        if (plugin->filename_extension) {
            char *buf = strdup(plugin->filename_extension);
            const char *token=strtok(buf, ",");
            if (token) do {
                extensiondict().insert(
                        ExtensionDict::value_type(std::string(token),plugin));
            } while ( (token = strtok(NULL, ",")) );
            free(buf);
        }
        return MOLFILE_SUCCESS;
    }

    struct _ {
        _() {
            MOLFILE_INIT_ALL;
            MOLFILE_REGISTER_ALL(NULL, register_cb);
        }
    } _static_initializer;

}

namespace desres { namespace molfile {


    const molfile_plugin_t *plugin_for_type(const char *type) {
        if (plugindict().count(type)>0) return plugindict()[type];
        return NULL;
    }
    
    
    const molfile_plugin_t *plugin_for_path(const char * path,
                                            const char * fallback ) {
        if (!path) return NULL;
        const char * ext = strrchr(path, '.');
        const molfile_plugin_t * plugin = NULL;
        if (ext) {
            ext += 1;
            if (extensiondict().count(ext)>0)
                plugin = extensiondict().find(ext)->second;
        }
        if (!plugin && fallback!=NULL && extensiondict().count(fallback)>0)
            plugin = extensiondict().find(fallback)->second;
        return plugin;
    }

    Frame::Frame(size_t natoms, bool with_velocities, bool double_precision)
    : m_natoms(natoms) {
        memset(&ts, 0, sizeof(ts));
        if (double_precision) {
            ts.dcoords = new double[3*natoms];
            memset(ts.dcoords, 0, 3*natoms*sizeof(double));
            if (with_velocities) {
                ts.dvelocities = new double[3*natoms];
                memset(ts.dvelocities, 0, 3*natoms*sizeof(double));
            }
        } else {
            ts.coords = new float[3*natoms];
            memset(ts.coords, 0, 3*natoms*sizeof(float));
            if (with_velocities) {
                ts.velocities = new float[3*natoms];
                memset(ts.velocities, 0, 3*natoms*sizeof(float));
            }
        }

        ts.total_energy = HUGE_VAL;
        ts.kinetic_energy = HUGE_VAL;
        ts.potential_energy = HUGE_VAL;
        ts.extended_energy = HUGE_VAL;
        ts.temperature = HUGE_VAL;
        ts.pressure = HUGE_VAL;
    }
    Frame::~Frame() {
        delete [] ts.coords;
        delete [] ts.velocities;
        delete [] ts.dcoords;
        delete [] ts.dvelocities;
    }

    Reader::~Reader() {
        if (plugin->close_file_read) plugin->close_file_read(handle);
    }

    Reader::Reader(const molfile_plugin_t *p, const char * _path,
                   bool double_precision) 
    : plugin(p), path(_path), handle(NULL), m_nframes(-1), m_natoms(0), 
      m_optflags(0), m_has_velocities(false), 
      m_double_precision(double_precision)
    {
            int natoms=0;
            handle = plugin->open_file_read(_path, plugin->name, &natoms);
            if (!handle) throw std::runtime_error("open_file_read failed");
            if (natoms<0) throw std::runtime_error("plugin could not determined #atoms");
            m_natoms = natoms;

            // count the number of frames
            bool supports_doubles = false;
            if (plugin->read_timestep_metadata) {
                molfile_timestep_metadata_t meta[1];
                memset(meta, 0, sizeof(meta));
                if (plugin->read_timestep_metadata(handle, meta)==MOLFILE_SUCCESS) {
                    m_nframes = meta->count;
                    m_has_velocities = meta->has_velocities;
                    supports_doubles = meta->supports_double_precision;
                }
            }
            /* plugin must explicitly declare its support for reading
             * of double precision positions and velocities */
            if (double_precision && !supports_doubles) {
                throw std::runtime_error(
                        std::string("plugin '") + p->name + "' does not support reading of double-precision positions"); 
            }

            // read structure information if available
            if (plugin->read_structure) {
                m_atoms.resize(natoms);
                if (plugin->read_structure(handle, &m_optflags, &m_atoms[0])
                        !=MOLFILE_SUCCESS) {
                    throw std::runtime_error("read_structure failed");
                }
                // read the bonds, if possible
                int *from=NULL, *to=NULL;
                float *order=NULL;
                int *bondtype=NULL;
                int nbondtypes;
                char **bondtypename;
                int nbonds=0;
                if (plugin->read_bonds && 
                        plugin->read_bonds(handle, &nbonds, &from, &to, &order,
                            &bondtype, &nbondtypes, &bondtypename) 
                        != MOLFILE_SUCCESS) { 
                    throw std::runtime_error("Failed reading bonds");
                }
                for (int i=0; i<nbonds; i++) {
                    ssize_t ai=from[i]-1;
                    ssize_t aj=to[i]-1;
                    if (ai<0 || aj<0 || ai>=natoms || aj>=natoms) {
                        throw std::runtime_error("Bond out of range");
                    }
                    m_bonds.push_back(bond_t(ai, aj, order ? order[i] : 1));
                }
            }
            if (plugin->read_volumetric_metadata) {
                int nsets;
                molfile_volumetric_t* meta=NULL;
                if (plugin->read_volumetric_metadata(handle, &nsets, &meta)
                        != MOLFILE_SUCCESS) {
                    throw std::runtime_error("Failed reading grid metadata");
                }
                m_grids.resize(nsets);
                if (nsets) {
                    memcpy(&m_grids[0], meta, nsets*sizeof(*meta));
                }
            }
        }

    Reader* Reader::reopen() const {
        return new Reader(plugin, path.c_str(), m_double_precision);
    }

    ssize_t Reader::read_times(ssize_t start, ssize_t count, double * times) const {
        if (!plugin->read_times) return -1;
        return plugin->read_times(handle, start, count, times);
    }

    namespace {
        struct Oracle {
            const Reader& reader;
            Oracle(const Reader& _reader) : reader(_reader) {}
            double operator[](ssize_t i) const {
                double T;
                if (reader.read_times(i,1,&T)!=1) {
                    throw std::runtime_error("Error reading time");
                }
                return T;
            }
        };
    }

#ifndef WIN32
    ssize_t Reader::Reader::at_time_near(double T) const {
        return ff::at_time_near<Oracle>(nframes(), *this, T);
    }
#endif
    ssize_t Reader::at_time_gt(double T) const {
        return ff::at_time_gt<Oracle>(nframes(), *this, T);
    }
    ssize_t Reader::at_time_ge(double T) const {
        return ff::at_time_ge<Oracle>(nframes(), *this, T);
    }
    ssize_t Reader::at_time_lt(double T) const {
        return ff::at_time_lt<Oracle>(nframes(), *this, T);
    }
    ssize_t Reader::at_time_le(double T) const {
        return ff::at_time_le<Oracle>(nframes(), *this, T);
    }

    Frame *Reader::frame(ssize_t index) const {
        if (index < 0) {
            if (m_nframes<0) {
                throw std::runtime_error(
                        "Cannot use negative index when number of frames is not known");
            }
            index += m_nframes;
        }
        Frame *result = new Frame(natoms(), has_velocities(), 
                                  double_precision());
        if (read_frame(index, *result) != MOLFILE_SUCCESS) {
            delete result;
            throw std::runtime_error("Reading frame failed");
        }
        return result;
    }

    int Reader::read_frame(ssize_t index, molfile_timestep_t* ts) const {
        if (!plugin->read_timestep2) 
            throw std::runtime_error("frame() not implemented for this plugin");
        return plugin->read_timestep2(handle, index, ts);
    }

    Frame *Reader::next() const {
        if (!plugin->read_next_timestep) return NULL;
        Frame *result = new Frame(natoms(), has_velocities(), double_precision());
        if (plugin->read_next_timestep(handle, natoms(), *result)!=MOLFILE_SUCCESS) {
            delete result;
            result=NULL;
        }
        return result;
    }

    void Reader::skip() const {
        if (plugin->read_next_timestep) {
            plugin->read_next_timestep(handle, natoms(), NULL);
        }
    }

    void Reader::read_grid(int n, float* data) const {
        if (n<0 || n>=(int)m_grids.size()) {
            throw std::runtime_error("Reader::read_grid - invalid grid");
        }
        if (!plugin->read_volumetric_data) {
            throw std::runtime_error("Reader::read_grid - plugin does not support reading grid data");
        }
        if (plugin->read_volumetric_data(handle, n, data, NULL)!=MOLFILE_SUCCESS) {
            throw std::runtime_error("Reader::read_grid - read failed");
        }
    }

    Writer::Writer(const molfile_plugin_t *p, const char * path, ssize_t natoms)
        : plugin(p), handle(NULL), m_natoms(natoms) {
            handle = plugin->open_file_write(path, plugin->name, natoms);
            if (!handle) throw std::runtime_error("Failed opening file for writing");
        }

    void Writer::close() {
        if (handle) plugin->close_file_write(handle);
        handle = NULL;
    }

    bool Writer::can_sync() const {
        return plugin->sync_file_write != NULL;
    }

    bool Writer::truncate(double t) {
        return plugin->truncate_file_write && 
               plugin->truncate_file_write(handle,t);
    }

    void Writer::sync() {
        if (plugin->sync_file_write) plugin->sync_file_write(handle);
    }

    Writer::~Writer() {
        close();
    }

    bool Writer::requires_atoms() const {
        return plugin->write_structure;
    }

    Writer& Writer::write_atoms(const std::vector<atom_t> &atoms,
                                const std::vector<bond_t> &bonds,
                                int optflags) {
        if (!handle) throw std::runtime_error("I/O on closed writer.");
        if (!requires_atoms()) return *this;
        
        std::vector<int> from, to;
        std::vector<float> order;
        if (plugin->write_bonds) {
            for (unsigned i=0; i<bonds.size(); i++) {
                from.push_back(bonds[i].from+1);
                to.push_back(bonds[i].to+1);
                order.push_back(bonds[i].order);
            }
            /* prevent illegal reference to element zero in the case 
             * of no bonds. */
            if (!bonds.size()) {
                from.resize(1);
                to.resize(1);
                order.resize(1);
            }
            if (plugin->write_bonds(handle, bonds.size(), 
                        &from[0], &to[0], &order[0], 
                        NULL, 0, NULL) != MOLFILE_SUCCESS) {
                throw std::runtime_error("Writing bonds failed");
            }
        }
        int rc;
        if (!atoms.size()) {
            /* provide dummy atoms */
            std::vector<atom_t> dummies(natoms());
            atom_t * ptr = natoms() ? &dummies[0] : NULL;
            memset(ptr, 0, natoms()*sizeof(*ptr));
            rc=plugin->write_structure(handle, 0, ptr);
        } else {
            rc=plugin->write_structure(handle, optflags, &atoms[0]);
        }
        if (rc != MOLFILE_SUCCESS) {
            throw std::runtime_error("Writing structure failed");
        }
        return *this;
    }

    bool Writer::writes_frames() const {
        return plugin->write_timestep;
    }

    Writer& Writer::write_frame(const Frame &frame) {
        if (!writes_frames())
            throw std::runtime_error("Plugin does not support writing frames");
        if ((ssize_t)frame.natoms() != natoms())
            throw std::runtime_error("Number of atoms in frame doesn't match file");
        if (!handle) throw std::runtime_error("I/O on closed writer.");
        if (plugin->write_timestep(handle, frame) != MOLFILE_SUCCESS)
            throw std::runtime_error("Writing frame failed");
        return *this;
    }

    Writer& Writer::write_grid(grid_t const& meta, const float* data) {
        if (!plugin->write_volumetric_data) {
            throw std::runtime_error("Plugin does not support writing grids");
        }
        /* FIXME: Why doesn't the C API take const pointers? */
        molfile_volumetric_t v(meta);
        float* d = const_cast<float *>(data);
        if (plugin->write_volumetric_data(handle, &v, d, NULL)) {
            throw std::runtime_error("Writing grid failed");
        }
        return *this;
    }

}}
