#include "../ff.hxx"

namespace {

    struct InPlaneWag : public Ffio {

        void apply( SystemPtr h,
                    const Json& blk,
                    const SiteMap& sitemap,
                    const VdwMap&, bool alchemical  ) const {

            TermTablePtr table = AddTable(h, "inplanewag_harm");
            ParamMap map(table->params(), blk);
            const Json& ai = blk.get("ffio_ai");
            const Json& aj = blk.get("ffio_aj");
            const Json& ak = blk.get("ffio_ak");
            const Json& al = blk.get("ffio_al");
            const Json& fn = blk.get("ffio_funct");
            IdList ids(4);
            int i,n = blk.get("__size__").as_int();
            for (i=0; i<n; i++) {
                std::string f = fn.elem(i).as_string();
                boost::to_lower(f);
                if (f!="harm") {
                    FFIO_ERROR("Expected ffio_funct='harm' in ffio_inplanewags; got " << f);
                }
                ids[0]=ai.elem(i).as_int();
                ids[1]=aj.elem(i).as_int();
                ids[2]=ak.elem(i).as_int();
                ids[3]=al.elem(i).as_int();
                sitemap.addUnrolledTerms( table, map.add(i), ids );
            }
        }
    };

    RegisterFfio<InPlaneWag> _("ffio_inplanewags");
}

