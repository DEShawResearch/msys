#ifndef MOLFILE_MOLFILE_HXX
#define MOLFILE_MOLFILE_HXX

#include "molfile_plugin.h"

#include <string>
#include <vector>
#include <cstdlib>
#include <cmath> // for HUGE_VAL

// This file contains C++ wrappers for the structure and coordinate reading
// parts of the VMD molfile_plugin API.  

namespace desres { namespace molfile {

    /* wrapper for molfile_timestep_t. */
    class Frame {
        molfile_timestep_t ts;
        const size_t m_natoms;

    public:
        Frame(size_t natoms, bool with_velocities=true,
                             bool double_precision=false);
        ~Frame();

        // act like a molfile_timestep_t so that we can pass it directly to
        // molfile_plugin API methods.
        operator       molfile_timestep_t*()       { return &ts; }
        operator const molfile_timestep_t*() const { return &ts; }

        size_t natoms() const { return m_natoms; }

        /* single-precision accessors when not using double precision */
        float *pos() { return ts.coords; }
        float *vel() { return ts.velocities; }
        const float *pos() const { return ts.coords; }
        const float *vel() const { return ts.velocities; }

        /* for consistency, single-precision with extended naming convention */
        float *fpos() { return pos(); }
        float *fvel() { return vel(); }
        const float *fpos() const { return pos(); }
        const float *fvel() const { return vel(); }

        /* double-precision accessors are required when the Frame is 
         * initialized with double_precision=true */
        double* dpos() { return ts.dcoords; }
        double* dvel() { return ts.dvelocities; }
        const double* dpos() const { return ts.dcoords; }
        const double* dvel() const { return ts.dvelocities; }

        double time() const { return ts.physical_time; }
        void setTime(double t) { ts.physical_time=t; }
        double *box() { return ts.unit_cell; }
        const double *box() const { return ts.unit_cell; }

        double* pressure_tensor() { return ts.pressure_tensor; }
        const double* pressure_tensor() const { return ts.pressure_tensor; }

        double* virial_tensor() { return ts.virial_tensor; }
        const double* virial_tensor() const { return ts.virial_tensor; }


#define MOLFILE_SCALAR_ACCESSOR(x) \
        bool has_##x () const { return ts. x != HUGE_VAL; } \
        double x () const { return ts. x ; } \
        void set_##x ( double v ) { ts. x = v; } 

        MOLFILE_SCALAR_ACCESSOR(total_energy)
        MOLFILE_SCALAR_ACCESSOR(kinetic_energy)
        MOLFILE_SCALAR_ACCESSOR(potential_energy)
        MOLFILE_SCALAR_ACCESSOR(extended_energy)
        MOLFILE_SCALAR_ACCESSOR(temperature)
        MOLFILE_SCALAR_ACCESSOR(pressure)

#undef MOLFILE_SCALAR_ACCESSOR
        // TODO: add select() and moveby()
    };

    typedef molfile_atom_t atom_t;

    struct bond_t {
        size_t from, to;
        float order;
        bond_t() {}
        bond_t(size_t f, size_t t, float o)
            : from(f), to(t), order(o) {}
    };

    typedef molfile_volumetric_t grid_t;

    class Reader {
        const molfile_plugin_t *plugin;
        std::string path;
        void *handle;
        ssize_t m_nframes;
        ssize_t m_natoms;
        std::vector<atom_t> m_atoms;
        std::vector<bond_t> m_bonds;
        std::vector<grid_t> m_grids;
        int m_optflags;
        bool m_has_velocities;
        bool m_double_precision;

    public:
        Reader(const molfile_plugin_t *p, const char * path,
               bool double_precision = false);
        ~Reader();

        /* create a new reader for the original file */
        Reader* reopen() const;

        const std::vector<atom_t>& atoms() const { return m_atoms; }
        const std::vector<bond_t>& bonds() const { return m_bonds; }
        const std::vector<grid_t>& grids() const { return m_grids; }

        /* read count times beginning at start.  Return number of times
         * read, or -1 if reading times isn't supported by the plugin. */
        ssize_t read_times(ssize_t start, ssize_t count, double * times) const;

        /* return index of corresponding frame, or -1 if no such frame could
         * be found. */
        ssize_t at_time_near(double T) const;
        ssize_t at_time_gt(double T) const;
        ssize_t at_time_ge(double T) const;
        ssize_t at_time_lt(double T) const;
        ssize_t at_time_le(double T) const;

        std::vector<atom_t>& atoms() { return m_atoms; }
        std::vector<bond_t>& bonds() { return m_bonds; }

        ssize_t nframes() const { return m_nframes; }
        ssize_t natoms() const { return m_natoms; }
        int optflags() const { return m_optflags; }
        bool has_velocities() const { return m_has_velocities; }
        bool double_precision() const { return m_double_precision; }
        Frame *frame(ssize_t index) const;

        int read_frame(ssize_t index, molfile_timestep_t* ts) const;
        Frame *next() const;
        void skip() const;

        /* Read n'th grid data into given space */
        void read_grid(int n, float* data) const;
    };

    class Writer {
        const molfile_plugin_t *plugin;
        void *handle;
        std::string m_path;
        ssize_t m_natoms;
    public:
        Writer(const molfile_plugin_t *p, const char * path, ssize_t natoms);
        ~Writer();

        bool can_sync() const;
        void sync();

        bool truncate(double after_time);

        void close();

        std::string path() const { return m_path; }

        bool requires_atoms() const;
        bool writes_frames() const;
        ssize_t natoms() const { return m_natoms; }

        Writer& write_atoms(const std::vector<atom_t> &atoms, 
                            const std::vector<bond_t> &bonds,
                            int optflags);
        Writer& write_frame(const Frame &frame);

        Writer& write_grid(grid_t const& meta, const float* data);
    };

    /*!
     * Return an instance of the plugin of the given type, or NULL if 
     * none is found.
     */
    const molfile_plugin_t * plugin_for_type(const char *type);
    
    /*!
     * Return an instance of the plugin for the given path.  If the type
     * could not be determined, use the provided fallback type.
     */
    const molfile_plugin_t *plugin_for_path(const char * path,
                                            const char * fallback = NULL );

}}

#endif
