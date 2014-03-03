#ifndef desres_msys_cealign_hxx
#define desres_msys_cealign_hxx

#include <vector>
#include <cmath>
#include <stdio.h>
#include <cassert>
#include <algorithm>

#include <pfx/rms.hxx>
#include <pfx/cell.hxx>

#include "types.hxx"

namespace desres { namespace msys {

    template <typename scalar>
    class CEAlign {
        typedef std::pair<unsigned, unsigned> AFP;
        typedef std::vector<AFP> Path;
        typedef std::vector<Path> PathList;

        const unsigned m;   /* window size */
        const unsigned G;   /* max gap between AFPs */
        const scalar D0;    /* threshold for single AFP */
        const scalar D1;    /* threshold for average pairwise AFPs */

        typedef std::vector<scalar> Vec;

        static
        Vec calc_distance_matrix(IdList const& inds,
                                 const scalar* pos) {

            unsigned n = inds.size();
            Vec dm(n*n);
            for (unsigned i=0; i<n; i++) {
                const scalar* pi = pos+3*inds[i];
                for (unsigned j=0; j<n; j++) {
                    const scalar* pj = pos+3*inds[j];
                    scalar d=0;
                    for (int k=0; k<3; k++) {
                        scalar delta = pi[k] - pj[k];
                        d += delta*delta;
                    }
                    d = std::sqrt(d);
                    dm[n*i+j] = d;
                }
            }
            return dm;
        }

#define DMA(_i, _j) (dmA.at(nA*(_i)+(_j)))
#define DMB(_i, _j) (dmB.at(nB*(_i)+(_j)))

        Vec calc_similarities(unsigned nA, Vec const& dmA, 
                              unsigned nB, Vec const& dmB) const {

            assert(nA*nA==dmA.size());
            assert(nB*nB==dmB.size());

            Vec S(nA*nB, D0);
            const unsigned sz = ((m-1)*(m-2))/2;

            for (unsigned i=0; i<nA; i++) {
                if (i+m>nA) continue;
                for (unsigned j=0; j<nB; j++) {
                    if (j+m>nB) continue;

                    scalar score=0;
                    /* diagonal entries in dm are of course 0; bidiagonal
                     * entries are all approximately equal, so dmA()-dmB()
                     * is close to zero so we can ignore those as well. */
                    for (unsigned k=0; k<m-2; k++) {
                        for (unsigned l=k+2; l<m; l++) {
                            score += std::fabs(DMA(i+k,i+l) - DMB(j+k,j+l));
                        }
                    }
                    S[i*nB+j] = score/sz;
                    //printf("%d %d %f\n",i,j,score/sz);
                }
            }
            return S;
        }

        PathList find_paths( unsigned nA, Vec const& dmA, 
                             unsigned nB, Vec const& dmB,
                             Vec const & S) const {

            PathList paths;
            unsigned best_path_length=1;

            for (unsigned iA=0; iA<nA; iA++) {
                /* stop if we can't match the longest path yet found */
                if (iA + m*(best_path_length-1) > nA) break;

                for (unsigned iB=0; iB<nB; iB++) {
                    /* stop if we can't match the longest path yet found */
                    if (iB+m*(best_path_length-1) > nB) break;

                    /* consider only plausible AFPs */
                    if (S[iA*nB+iB] >= D0) continue;

                    /* starting a new path */
                    Path path(1, AFP(iA,iB));

                    for (;;) {
                        /* consider all possible gaps */
                        AFP afp;
                        scalar gap_best_score = D1;     /* worst possible */
                        for (unsigned g=0; g<2*G+1; g++) {
                            unsigned jA = m+path.back().first;
                            unsigned jB = m+path.back().second;
                            if (g%2) jA += (g+1)/2;
                            else     jB += (g+1)/2;
                            if (jA+m >= nA || jB+m >= nB) continue;
                            if (S[jA*nB+jB] >= D0) continue;

                            /* compute score for proposed extension */
                            scalar score=0;
                            for (unsigned i=0, n=path.size(); i<n; i++) {
                                unsigned ia = path[i].first;
                                unsigned ib = path[i].second;
                                scalar s1 = fabs( DMA(ia,     jA) - 
                                                  DMB(ib,     jB) );
                                scalar sm = fabs( DMA(ia+m-1, jA+m-1) -
                                                  DMB(ib+m-1, jB+m-1) );
                                scalar sk=0;
                                for (unsigned k=1; k<m-1; k++) {
                                    sk +=   fabs( DMA(ia+k,   jA+(m-1)-k) -
                                                  DMB(ib+k,   jB+(m-1)-k) );
                                }
                                score += s1 + sm + sk;
                            }
                            score /= m*path.size();
                            if (score < gap_best_score) {
                                gap_best_score = score;
                                afp = AFP(jA, jB);
                            }
                        }
                        if (gap_best_score<D1) {
                            path.push_back(afp);
                            continue;
                        }
                        break;
                    }

                    /* keep path if same length or longer than current best */
                    if (path.size() >= best_path_length) {
                        best_path_length = path.size();
                        paths.push_back(path);
                    }
                }   /* loop over iB */
            }   /* loop over iA */

            return paths;
        }
        static
        double find_best(PathList const& paths, unsigned m,
                         const unsigned* A, const scalar* Apos,
                         const unsigned* B, const scalar* Bpos,
                         scalar* mat_3x3, scalar* T1, scalar* T2) {
            if (paths.empty()) return 0;
            double best_rmsd=HUGE_VAL;
            scalar mat[9];
            unsigned sz=paths.back().size();        /* number of AFPs */
            /* consider only the longest paths */
            for (unsigned i=paths.size(); i--;) {
                Path const& path = paths[i];
                if (path.size()<sz) break;
                std::vector<scalar> ref, pos;
                ref.reserve(sz*3*m);
                pos.reserve(sz*3*m);
                for (unsigned j=0; j<sz; j++) {
                    AFP const& afp = path[j];
                    unsigned a0 = afp.first;
                    unsigned b0 = afp.second;
                    for (unsigned k=0; k<m; k++) {
                        const scalar* pa = Apos + 3*A[a0+k];
                        const scalar* pb = Bpos + 3*B[b0+k];
                        ref.insert(ref.end(), pa, pa+3);
                        pos.insert(pos.end(), pb, pb+3);
                    }
                }
                scalar tmp1[3], tmp2[3];
                pfx::compute_center(m*sz, NULL, &ref[0], tmp1);
                pfx::apply_shift(m*sz, &ref[0], -tmp1[0], -tmp1[1], -tmp1[2]);
                pfx::compute_center(m*sz, NULL, &pos[0], tmp2);
                pfx::apply_shift(m*sz, &pos[0], -tmp2[0], -tmp2[1], -tmp2[2]);
                double rmsd = pfx::compute_alignment(m*sz, NULL,
                                                &ref[0], &pos[0], mat);
                if (rmsd<best_rmsd) {
                    best_rmsd=rmsd;
                    memcpy(mat_3x3, mat, sizeof(mat));
                    memcpy(T1, tmp1, sizeof(tmp1));
                    memcpy(T2, tmp2, sizeof(tmp2));
                }
            }
            return best_rmsd;
        }

    public:
        /* Construct a CEAlign engine.
         * The standard defaults are:
         *      window_size=8
         *      max_gap=30
         *      local_threshold=3.0
         *      pair_threshold=4.0
         * but all values must be supplied here.
         */
        CEAlign(unsigned window_size, unsigned max_gap,
                scalar local_threshold, scalar pair_threshold)
        : m(window_size) 
        , G(max_gap)
        , D0(local_threshold)
        , D1(pair_threshold)
        {}

        unsigned window_size() const { return m; }

        /* Construct a CEAlign engine with the published defaults */
        static CEAlign<scalar> WithDefaults() {
            return CEAlign(8, 30, scalar(3.0), scalar(4.0));
        }

        typedef std::pair<IdList, IdList> Map;
        typedef std::vector<Map> MapList;

        Map path_to_map(Path const& path) const {
            Map map;
            for (Path::const_iterator i=path.begin(), e=path.end(); i!=e; ++i) {
                unsigned a0 = i->first;
                unsigned b0 = i->second;
                for (unsigned k=0; k<m; k++) {
                    map.first.push_back(a0+k);
                    map.second.push_back(b0+k);
                }
            }
            return map;
        }

        /* Find good mappings of atoms B onto atoms in A. */
        MapList compute(IdList const& A, const scalar* Apos,
                        IdList const& B, const scalar* Bpos) const {

            /* precalculate pairwise distances */
            Vec dmA = calc_distance_matrix(A, Apos);
            Vec dmB = calc_distance_matrix(B, Bpos);

            /* precalculate AFP similarities using metric (ii) */
            Vec S = calc_similarities(A.size(), dmA, B.size(), dmB);

            /* find the paths */
            PathList paths = find_paths(A.size(), dmA, B.size(), dmB, S);
            if (paths.empty()) return MapList();

            /* keep only the longest paths */
            size_t longest = paths.back().size();
            PathList::iterator it = paths.begin(), end=paths.end();
            while (it->size() < longest) ++it;

            /* convert paths into full mappings */
            MapList maps;
            for (; it!=end; ++it) maps.push_back(path_to_map(*it));
            return maps;
        }
    };

}}


#endif
