#include "../ff.hxx"

namespace {

    static const char* angle_cols[] = {"ffio_sigma", "ffio_theta0", "ffio_fc"};
    static const char* improper_cols[] = {"ffio_sigma", "ffio_phi0", "ffio_fc"};
    static const char* posre_cols[] = {"ffio_fc", "ffio_sigma"};

    struct Angle : public Ffio {

        void apply( SystemPtr h,
                    const Json& blk,
                    const SiteMap& sitemap,
                    const VdwMap& ) const {

            TermTablePtr table = AddTable(h,"angle_fbhw");
            ParamMap map(table->params(), blk, 3, angle_cols);

            const Json& ai = blk.get("ffio_ai");
            const Json& aj = blk.get("ffio_aj");
            const Json& ak = blk.get("ffio_ak");
            int i,n = blk.get("__size__").as_int();
            IdList ids(3);
            for (i=0; i<n; i++) {
                ids[0]=ai.elem(i).as_int()-1;
                ids[1]=aj.elem(i).as_int()-1;
                ids[2]=ak.elem(i).as_int()-1;
                table->addTerm(ids, map.add(i));
            }
        }
    };

    struct Improper : public Ffio {

        void apply( SystemPtr h,
                    const Json& blk,
                    const SiteMap& sitemap,
                    const VdwMap& ) const {

            TermTablePtr table = AddTable(h,"improper_fbhw");
            ParamMap map(table->params(), blk, 3, improper_cols);

            const Json& ai = blk.get("ffio_ai");
            const Json& aj = blk.get("ffio_aj");
            const Json& ak = blk.get("ffio_ak");
            const Json& al = blk.get("ffio_al");
            int i,n = blk.get("__size__").as_int();
            IdList ids(4);
            for (i=0; i<n; i++) {
                ids[0]=ai.elem(i).as_int()-1;
                ids[1]=aj.elem(i).as_int()-1;
                ids[2]=ak.elem(i).as_int()-1;
                ids[3]=al.elem(i).as_int()-1;
                table->addTerm(ids, map.add(i));
            }
        }
    };


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

    RegisterFfio<Posre>     _1("ffio_posre_fbhw");
    RegisterFfio<Angle>     _2("ffio_angle_fbhw");
    RegisterFfio<Improper>  _3("ffio_improper_fbhw");
}

