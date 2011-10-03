#ifndef desres_msys_mae_sitemap_hxx
#define desres_msys_mae_sitemap_hxx

#include "../system.hxx"
#include "../term_table.hxx"
#include <fastjson/fastjson.hxx>

namespace desres { namespace msys { namespace mae {

    class SiteMap {
        IdList  _atoms;
        Id      _nsites;
        std::vector<Id> _s2p;

    public:
        SiteMap( SystemPtr h, const fastjson::Json& sites, const IdList& atoms, 
                 int natoms, int npseudos );

        Id nsites() const { return _nsites; }

        const IdList& atoms() const { return _atoms; }
        Id site(Id id) const {
            return _atoms.at(_s2p[id-1]);
        }

        Id addTerm(TermTablePtr table, 
                     Id p,
                     const IdList& ids) const;

        void addUnrolledTerms( TermTablePtr nb, 
                               Id param, 
                               const IdList& ids,
                               bool constrained, 
                               Id paramB ) const;

        void addUnrolledTerms( TermTablePtr nb, 
                               Id param, 
                               const IdList& ids,
                               bool constrained = false ) const {
            return addUnrolledTerms( nb, param, ids, constrained, BadId );
        }

        void addUnrolledTerms( TermTablePtr nb, 
                               Id param, 
                               const IdList& ids,
                               Id paramB ) const {
            return addUnrolledTerms( nb, param, ids, false, paramB );
        }

        void addUnrolledPseudoBonds(SystemPtr h, Id parent, Id pseudo) const;
    };

}}}

#endif
