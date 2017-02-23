#include "../ff.hxx"

namespace {

    /* return true if matches harm(_constrained)+ */
    bool harm_constrained(const std::string& s) {
        return s.substr(0,5)=="harm_" && (s.size()-s.rfind("_constrained"))==12;
    }

    template <bool alchemical>
    struct Bonds : public Ffio {

#ifdef DESMOND_USE_SCHRODINGER_MMSHARE
        typedef std::vector<const char *> StrList;
        static StrList make_columns(int n) {
            static const char * regcols[] = {
                "ffio_c1", "ffio_c2", "ffio_c3", "ffio_c4"
            };
            static const char * alccols[] = {
                "ffio_c1B", "ffio_c2B", "ffio_c3B", "ffio_c4B"
            };
            StrList cols(regcols,regcols+n);
            if (alchemical) cols.insert(cols.end(), alccols, alccols+n);
            return cols;
        }
#endif

        void apply( SystemPtr h,
                    const Json& blk,
                    const SiteMap& sitemap,
                    const VdwMap& ) const {

#ifdef DESMOND_USE_SCHRODINGER_MMSHARE
            std::string prefix;
            if (alchemical) prefix = "alchemical_";

            TermTablePtr bonds, softbonds;
            ParamMapPtr hmap, smap;

            const Json& sc = blk.get("ffio_schedule");
            
            TermTablePtr table;
#else
            std::string name;
            if (alchemical) name += "alchemical_";
            name += "stretch_harm";
            TermTablePtr table = AddTable(h, name);
            ParamMap map(table->params(), blk, alchemical);
#endif
            const Json& ai = blk.get("ffio_ai");
            const Json& aj = blk.get("ffio_aj");
            const Json& fn = blk.get("ffio_funct");

            IdList ids(2);

            int i,n = blk.get("__size__").as_int();
            for (i=0; i<n; i++) {
                std::string f = fn.elem(i).as_string();
                to_lower(f);
#ifdef DESMOND_USE_SCHRODINGER_MMSHARE
                const char* s = sc.valid() ? sc.elem(i).as_string() : 0;

                Id A=BadId;

                bool constrained = harm_constrained(f);
                if (f=="harm" || constrained) {
                    if (!bonds) {
                        bonds = AddTable(h, prefix+"stretch_harm");
                        StrList cols = make_columns(2);
                        hmap.reset(
                            new ParamMap(bonds->params(), blk,
                                cols.size(), &cols[0]));
                    }
                    table = bonds;
                    A = hmap->add(i);
                } else if (f=="harm_soft") {
                    if (!softbonds) {
                        softbonds = AddTable(h, prefix+"softstretch_harm");
                        StrList cols = make_columns(3);
                        smap.reset(
                            new ParamMap(softbonds->params(), blk,
                                cols.size(), &cols[0]));
                    }
                    table = softbonds;
                    A = smap->add(i);
#else
                bool constrained=false;
                if (f=="harm") {
                } else if (harm_constrained(f)) {
                    constrained=true;
#endif
                } else {
                    FFIO_ERROR("Unsupported ffio_funct in ffio_bonds: " << f);
                }
#ifndef DESMOND_USE_SCHRODINGER_MMSHARE
                Id A = map.add(i);
#endif
                ids[0] = ai.elem(i).as_int();
                ids[1] = aj.elem(i).as_int();
                /* we can't really constrain alchemical terms... */
                if (alchemical && constrained) {
                    constrained = false;
                }
#ifdef DESMOND_USE_SCHRODINGER_MMSHARE
                sitemap.addUnrolledTerms( table, A, ids, constrained, s );
#else
                sitemap.addUnrolledTerms( table, A, ids, constrained );
#endif
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

