#include "../ff.hxx"

namespace {

    bool harm_constrained(const std::string& s) {
        return s.substr(0,5)=="harm_" && (s.size()-s.rfind("_constrained"))==12;
    }

    struct Angles : public Ffio {

        void apply( SystemPtr h,
                    const Json& blk,
                    const SiteMap& sitemap,
                    const VdwMap&, bool alchemical ) const {

            TermTablePtr bonds = AddTable(h,"stretch_harm");
            TermTablePtr angles= AddTable(h,"angle_harm");
            ParamMap bmap(bonds->params(), blk);
            ParamMap amap(angles->params(), blk);

            const Json& ai = blk.get("ffio_ai");
            const Json& aj = blk.get("ffio_aj");
            const Json& ak = blk.get("ffio_ak");
            const Json& fn = blk.get("ffio_funct");

            int i,n = blk.get("__size__").as_int();
            for (i=0; i<n; i++) {
                std::string f = fn.elem(i).as_string();
                boost::to_lower(f);
                bool constrained=false;
                TermTablePtr table;
                Id A=BadId, B=BadId;
                IdList ids;
                if (f=="ub") {
                    table=bonds;
                    A=bmap.add(i);
                    ids.push_back(ai.elem(i).as_int());
                    ids.push_back(ak.elem(i).as_int());
                } else if (f=="harm" || harm_constrained(f)) {
                    table=angles;
                    A=amap.add(i);
                    ids.push_back(ai.elem(i).as_int());
                    ids.push_back(aj.elem(i).as_int());
                    ids.push_back(ak.elem(i).as_int());
                    if (harm_constrained(f)) constrained=true;
                    if (alchemical) B = amap.addAlchemical(i);
                    /* we can't really constrain alchemical terms... */
                    if (alchemical && constrained) {
                        constrained = false;
                    }
                } else {
                    FFIO_ERROR("Unsupported ffio_funct in ffio_angles: " << f);
                }
                sitemap.addUnrolledTerms( table, A, ids, constrained, B );
            }
        }
    };

    RegisterFfio<Angles> _("ffio_angles");
}

