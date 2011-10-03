#include "../maeatoms.hxx"

namespace {

    void write_table( const Json& ff,
                      int cmapid,
                      SystemPtr h ) {

        char buf[32];
        sprintf(buf, "ffio_cmap%d", cmapid);
        const Json& blk = ff.get(buf);

        ParamTablePtr d(ParamTable::create());
        d->addProp("phi", FloatType);
        d->addProp("psi", FloatType);
        d->addProp("energy", FloatType);

        const Json& ai = blk.get("ffio_ai");
        const Json& aj = blk.get("ffio_aj");
        const Json& c1 = blk.get("ffio_c1");
        int i,n = blk.get("__size__").as_int();
        for (i=0; i<n; i++) {
            Id row = d->addParam();
            d->value(row,0)=ai.elem(i).as_float();
            d->value(row,1)=aj.elem(i).as_float();
            d->value(row,2)=c1.elem(i).as_float();
        }
        const char * name = buf+5;
        h->addExtra(name, d);
    }

    struct Cmap : public Ffio {

        virtual bool wants_all() const { return true; }

        void apply( SystemPtr h,
                    const Json& ff,
                    const SiteMap& sitemap,
                    const VdwMap&, bool alchemical  ) const {

            const Json& blk = ff.get("ffio_torsion_torsion");
            const Json& c1 = blk.get("ffio_c1");
            if (!c1) throw std::runtime_error(
                    "ffio_torsion_torsion missing ffio_c1");

            TermTablePtr table = AddTable(h,"torsiontorsion_cmap");
            ParamTablePtr params = table->paramTable();
            MaeAtoms atoms(blk);
            std::map<int,Id> pmap;

            int i,n = blk.get("__size__").as_int();
            for (i=0; i<n; i++) {
                int cmapid = c1.elem(i).as_int();
                if (!pmap.count(cmapid)) {
                    Id param = params->addParam();
                    char val[32];
                    sprintf(val, "cmap%d", cmapid);
                    params->value(param,0)=val;
                    pmap[cmapid]=param;
                    write_table(ff, cmapid, h);
                }
                sitemap.addUnrolledTerms(table, pmap[cmapid], atoms.ids(i));
            }
        }
    };

    struct Dummy : public Ffio {

        void apply( SystemPtr h,
                    const Json& blk,
                    const SiteMap& sitemap,
                    const VdwMap&, bool   ) const {
        }
    };

    RegisterFfio<Cmap> _("ffio_torsion_torsion");
    RegisterFfio<Dummy> _1("ffio_cmap1");
    RegisterFfio<Dummy> _2("ffio_cmap2");
    RegisterFfio<Dummy> _3("ffio_cmap3");
    RegisterFfio<Dummy> _4("ffio_cmap4");
    RegisterFfio<Dummy> _5("ffio_cmap5");
    RegisterFfio<Dummy> _6("ffio_cmap6");
}

