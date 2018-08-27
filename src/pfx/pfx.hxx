#ifndef desres_pfx_hxx
#define desres_pfx_hxx

#include "graph.hxx"
#include "cell.hxx"
#include "glue.hxx"
#include "rms.hxx"

#include "../system.hxx"

namespace desres { namespace msys { namespace pfx {

    struct Pfx {

        // topologically sorted bonds.  
        std::vector<unsigned> _bonds;

        // atoms sorted by connected component
        std::vector<unsigned> _comps;

        // size of each connected component
        std::vector<unsigned> _sizes;

        // glue elements
        std::vector<Glue>     _glue;

        // atoms involved in alignment
        std::vector<unsigned> _align;

        // sorted, unique align atoms for quick lookup
        std::vector<unsigned> _align_lookup;
        bool is_aligned(unsigned id) const {
            return std::binary_search(
                    _align_lookup.begin(),
                    _align_lookup.end(),
                    id);
        }

        // reference coordinates for alignment
        std::vector<double>   _aref;

        // weights for alignment.
        std::vector<double>   _weights;

        // center of reference structure
        double _cpos[3];

        // do periodic wrapping of bonds
        template <typename scalar>
        void fix_bonds(const scalar* cell, const scalar* proj,
                       scalar* pos) const {
            for (unsigned i=0,n = _bonds.size(); i<n; i+=2) {
                scalar* pi = pos+3*_bonds[i  ];
                scalar* pj = pos+3*_bonds[i+1];
                scalar d[3] = {pj[0]-pi[0], pj[1]-pi[1], pj[2]-pi[2]};
                wrap_vector(cell, proj, d);
                pj[0] += d[0];
                pj[1] += d[1];
                pj[2] += d[2];
            }
        }

        // do periodic wrapping of fragments
        template <typename scalar>
        void wrap_frags(const scalar* cell, const scalar* proj,
                        scalar* pos) const {
            if (size()==0) return;

            const unsigned *atom = &_comps[0];
            for (unsigned i=0, n=_sizes.size(); i<n; i++) {
                const unsigned sz = _sizes[i];
                // compute center of fragment
                scalar s = 1/(scalar)sz;
                scalar c[3]={0,0,0};
                bool skip = false;
                for (unsigned j=0; j<sz; j++) {
                    auto id = atom[j];
                    if (is_aligned(id)) {
                        skip = true;
                        break;
                    }
                    const scalar* p = pos+3*id;
                    c[0] += p[0];
                    c[1] += p[1];
                    c[2] += p[2];
                }
                if (!skip) {
                    c[0] *= s;
                    c[1] *= s;
                    c[2] *= s;
                    
                    // compute shift for wrapping
                    wrap_vector(cell, proj, c);

                    // apply shift to fragment
                    for (unsigned j=0; j<sz; j++) {
                        scalar* p = pos+3*atom[j];
                        p[0] += c[0];
                        p[1] += c[1];
                        p[2] += c[2];
                    }
                }
                // advance to next fragment
                atom += sz;
            }
        }

    public:
        Pfx(Graph const& g, bool fix_bonds=false) {

            if (fix_bonds) {
                g.copy_bonds(std::back_inserter(_bonds));
            }
            g.copy_components(std::back_inserter(_comps),
                              std::back_inserter(_sizes));
        }

        // number of atoms as defined by the initial graph
        unsigned size() const { return _comps.size(); }

        // Specify atoms to be wrapped as a single fragment.
        void glue(unsigned n, const unsigned* atoms) {
            _glue.push_back(Glue(_comps, _sizes, n, atoms));
        }

        // Specify that atoms should be oriented using the given reference 
        // coordinates, and centered on those coordinates. These atoms will 
        // be automatically added as an aggregate as well.  
        //
        // If coords is NULL, then the pipeline will center the specified atoms
        // at the origin.
        template <typename scalar>
        void align(unsigned n, const unsigned* atoms,
                               const scalar* coords,
                               const scalar* weights=nullptr) {

            _align.resize(n);
            std::copy(atoms, atoms+n, _align.begin());
            _align_lookup = _align;
            sort_unique(_align_lookup);
            if (weights) {
                _weights.resize(n);
                std::copy(weights, weights+n, _weights.begin());
            }
            const double* wts = _weights.empty() ? nullptr : &_weights[0];

            if (coords) {
                compute_center(n, nullptr, coords, _cpos, wts);
                _aref.resize(3*n);
                std::copy(coords, coords+3*n, _aref.begin());
                apply_shift(n, &_aref[0], -_cpos[0], -_cpos[1], -_cpos[2]);
            } else {
                _aref.clear();
            }
            glue(n, atoms);
        }

        // Perform previously specified wrapping and alignment operations
        // on the provided data.  pos must be of size Nx3 where N==size().
        // cell, if non-NULL, must be of size 3x3 and hold unit cell vectors
        // in C-major rows.  vel, if non-NULL, must be of size Nx3.
        template <typename scalar, typename cell_scalar>
        void apply(scalar* pos, cell_scalar* cell, scalar* vel) const {
            scalar box[9], inv[9];
            if (!pos) return;
            const double* wts = _weights.empty() ? nullptr : &_weights[0];
            if (cell) {
                trans_3x3(box, cell);
                if (!inverse_3x3(inv, box)) {
                    memset(inv, 0, sizeof(inv));
                }
                std::copy(cell, cell+9, box);

                // fix bonds
                fix_bonds(box, inv, pos);
                // apply glue
                for (unsigned i=0, n=_glue.size(); i<n; i++) {
                    _glue[i].join(box, inv, pos);
                }
            }
            // align or center
            if (!_align.empty()) {
                scalar center[3];
                compute_center(_align.size(), &_align[0], pos, center, wts);
                apply_shift(size(), pos, -center[0], -center[1], -center[2]);
                if (!_aref.empty()) {
                    // alignment
                    scalar mat[9], matT[9];
                    compute_alignment(_align.size(), &_align[0], &_aref[0],
                                      pos, mat, wts);
                    apply_rotation(size(), pos, mat);
                    if (vel) {
                        apply_rotation(size(), vel, mat);
                    }
                    if (cell) {
                        apply_rotation(3, box, mat);
                        trans_3x3(matT, mat);
                        // compute inv <- inv @ mat'
                        apply_rotation(3, matT, inv);
                        std::copy(matT, matT+9, inv);
                    }
                }
            }
            if (cell) {
                // wrap frags
                wrap_frags(box, inv, pos);

                // copy back to input
                std::copy(box, box+9, cell);

            }
            // final shift when doing alignment
            if (!_aref.empty()) {
                apply_shift<scalar>(size(), pos, _cpos[0], _cpos[1], _cpos[2]);
            }
        }

        // Compute rmsd with reference coordinates.  If none have been given
        // with pfx_align, return -1.
        template <typename scalar>
        double rmsd(const scalar* pos) const {
            if (_aref.empty()) return -1;
            const double* wts = _weights.empty() ? nullptr : &_weights[0];
            return compute_rmsd(_align.size(), &_align[0], 
                                &_aref[0], _cpos, pos, wts);
        }
    };

}}}

#endif
