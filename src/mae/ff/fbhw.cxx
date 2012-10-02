#include "../ff.hxx"

namespace {

    static const char* posre_cols[] = { "ffio_fc", "ffio_sigma" };

    struct Posre : public Ffio {

        void apply( SystemPtr h,
                    const Json& blk,
                    const SiteMap& sitemap,
                    const VdwMap& ) const {

            TermTablePtr table = AddTable(h,"posre_fbhw");
            ParamMap map(table->params(), blk, 2, posre_cols);

            const Json& ai = blk.get("ffio_ai");
            const Json& x0 = blk.get("ffio_x0");
            const Json& y0 = blk.get("ffio_y0");
            const Json& z0 = blk.get("ffio_z0");

            int i,n = blk.get("__size__").as_int();
            for (i=0; i<n; i++) {
                IdList ids(1,ai.elem(i).as_int()-1);
                ids[0] = sitemap.atoms().at(ids[0]);
                Id A = map.add(i);
                Id term = table->addTerm(ids, A);
                table->termPropValue(term,0)=x0.elem(i).as_float();
                table->termPropValue(term,1)=y0.elem(i).as_float();
                table->termPropValue(term,2)=z0.elem(i).as_float();
            }
        }
    };

    RegisterFfio<Posre> _("ffio_posre_fbhw");
}

