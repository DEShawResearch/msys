#include "../ff.hxx"

namespace {

    bool harm_constrained(const std::string& s) {
        return s.substr(0,5)=="harm_" && (s.size()-s.rfind("_constrained"))==12;
    }

    template <bool alchemical>
    struct Angles : public Ffio {

        void apply( SystemPtr h,
                    const Json& blk,
                    const SiteMap& sitemap,
                    const VdwMap& ) const {

            std::string bname;
            if (alchemical) bname = "alchemical_";
            bname += "stretch_harm";
            TermTablePtr bonds = AddTable(h,bname);

            std::string aname;
            if (alchemical) aname = "alchemical_";
            aname += "angle_harm";
            TermTablePtr angles= AddTable(h,aname);

            ParamMap bmap(bonds->params(), blk, alchemical);
            ParamMap amap(angles->params(), blk, alchemical);

            const Json& ai = blk.get("ffio_ai");
            const Json& aj = blk.get("ffio_aj");
            const Json& ak = blk.get("ffio_ak");
            const Json& fn = blk.get("ffio_funct");

            int i,n = blk.get("__size__").as_int();
            for (i=0; i<n; i++) {
                std::string f = fn.elem(i).as_string();
                to_lower(f);
                bool constrained=false;
                TermTablePtr table;
                Id A=BadId;
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
                    /* we can't really constrain alchemical terms... */
                    if (alchemical && constrained) {
                        constrained = false;
                    }
                } else {
                    FFIO_ERROR("Unsupported ffio_funct in ffio_angles: " << f);
                }
                sitemap.addUnrolledTerms( table, A, ids, constrained );
            }
        }
    };

    RegisterFfio<Angles<false> > _1("ffio_angles");
    RegisterFfio<Angles<true> > _2("ffio_angles_alchemical");
}

