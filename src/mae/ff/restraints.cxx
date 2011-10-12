#include "../ff.hxx"

namespace {

    struct Restraints : public Ffio {

        void apply( SystemPtr h,
                    const Json& blk,
                    const SiteMap& sitemap,
                    const VdwMap&, bool alchemical  ) const {

            TermTablePtr table = AddTable(h,"posre_harm");
            ParamMap map(table->params(), blk);

            const Json& ai = blk.get("ffio_ai");
            const Json& t1 = blk.get("ffio_t1");
            const Json& t2 = blk.get("ffio_t2");
            const Json& t3 = blk.get("ffio_t3");
            const Json& fn = blk.get("ffio_funct");

            int i,n = blk.get("__size__").as_int();
            for (i=0; i<n; i++) {
                std::string f = fn.elem(i).as_string("");
                boost::to_lower(f);
                if (f!="harm") {
                    FFIO_ERROR("Unsupported ffio_funct in ffio_restraints: " 
                            << f);
                }
                IdList ids(1,ai.elem(i).as_int()-1);
                ids[0] = sitemap.atoms().at(ids[0]);
                Id A = map.add(i);
                Id term = table->addTerm(ids, A);
                table->termPropValue(term,0)=t1.elem(i).as_float();
                table->termPropValue(term,1)=t2.elem(i).as_float();
                table->termPropValue(term,2)=t3.elem(i).as_float();
            }
        }
    };

    RegisterFfio<Restraints> _("ffio_restraints");
}

