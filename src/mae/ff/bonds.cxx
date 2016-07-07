#include "../ff.hxx"

namespace {

    /* return true if matches harm(_constrained)+ */
    bool harm_constrained(const std::string& s) {
        return s.substr(0,5)=="harm_" && (s.size()-s.rfind("_constrained"))==12;
    }

    template <bool alchemical>
    struct Bonds : public Ffio {

        void apply( SystemPtr h,
                    const Json& blk,
                    const SiteMap& sitemap,
                    const VdwMap& ) const {

            std::string name;
            if (alchemical) name += "alchemical_";
            name += "stretch_harm";
            TermTablePtr table = AddTable(h, name);
            ParamMap map(table->params(), blk, alchemical);

            const Json& ai = blk.get("ffio_ai");
            const Json& aj = blk.get("ffio_aj");
            const Json& fn = blk.get("ffio_funct");

            IdList ids(2);

            int i,n = blk.get("__size__").as_int();
            for (i=0; i<n; i++) {
                std::string f = fn.elem(i).as_string();
                to_lower(f);
                bool constrained=false;
                if (f=="harm") {
                } else if (harm_constrained(f)) {
                    constrained=true;
                } else {
                    FFIO_ERROR("Unsupported ffio_funct in ffio_bonds: " << f);
                }

                Id A = map.add(i);
                ids[0] = ai.elem(i).as_int();
                ids[1] = aj.elem(i).as_int();
                /* we can't really constrain alchemical terms... */
                if (alchemical && constrained) {
                    constrained = false;
                }
                sitemap.addUnrolledTerms( table, A, ids, constrained );
            }
        }
    };

    template <bool alchemical>
    struct Morse : public Ffio {

        void apply( SystemPtr h,
                    const Json& blk,
                    const SiteMap& sitemap,
                    const VdwMap& ) const {

            std::string name;
            if (alchemical) name += "alchemical_";
            name += "stretch_morse";
            TermTablePtr table = AddTable(h, name);
            ParamMap map(table->params(), blk, alchemical);

            IdList ids(2);

            const Json& ai = blk.get("ffio_ai");
            const Json& aj = blk.get("ffio_aj");

            int i,n = blk.get("__size__").as_int();
            for (i=0; i<n; i++) {
                Id A = map.add(i);
                ids[0] = ai.elem(i).as_int();
                ids[1] = aj.elem(i).as_int();
                sitemap.addUnrolledTerms( table, A, ids );
            }
        }
    };

    RegisterFfio<Bonds<false> > _1("ffio_bonds");
    RegisterFfio<Bonds<true> > _2("ffio_bonds_alchemical");
    RegisterFfio<Morse<false> > _3("ffio_morsebonds");
    RegisterFfio<Morse<true> > _4("ffio_morsebonds_alchemical");
}

