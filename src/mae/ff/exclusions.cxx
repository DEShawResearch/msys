#include "../ff.hxx"

namespace {

    struct Exclusions : public Ffio {
        void apply( SystemPtr h,
        const Json& blk,
        const SiteMap& sitemap,
        const VdwMap&, bool alchemical  ) const {

            TermTablePtr table = AddTable(h, "exclusion");
            IdList ids(2);
            const Json& ai = blk.get("ffio_ai");
            const Json& aj = blk.get("ffio_aj");
            int i,n = blk.get("__size__").as_int();
            for (i=0; i<n; i++) {
                ids[0]=ai.elem(i).as_int();
                ids[1]=aj.elem(i).as_int();
                sitemap.addUnrolledTerms( table, BadId, ids );
            }
        }
    };

    RegisterFfio<Exclusions> _("ffio_exclusions");
}

