#include "../ff.hxx"

namespace {

    static const char* angle_cols[] = {"ffio_sigma", "ffio_theta0", "ffio_fc"};
    static const char* improper_cols[] = {"ffio_sigma", "ffio_phi0", "ffio_fc"};
    static const char* posre_cols[] = {"ffio_fc", "ffio_sigma"};

    static IdList parse_ids(std::string const& s) {
        std::stringstream ss(s);
        IdList ids;
        Id id;
        while (ss >> id) ids.push_back(id-1);
        return ids;
    }

    struct Stretch : public Ffio {
        void apply( SystemPtr h,
                    const Json& blk,
                    const SiteMap& sitemap,
                    const VdwMap& ) const {

            static const char* stretch_cols[] = {
                "ffio_lower", "ffio_upper", "ffio_sigma", "ffio_beta", "ffio_fc"
            };

            Int max_group=0;
            TermTablePtr table = h->addTable("stretch_fbhw", 1);
            table->category = BOND;
            table->addTermProp("group", IntType);
            for (Id i=0; i<table->termCount(); i++) {
                max_group=std::min(
                        max_group, table->termPropValue(i,"group").asInt());
            }
            ParamTablePtr params = table->params();
            for (int i=0; i<5; i++) {
                params->addProp(stretch_cols[i]+5, FloatType);
            }
            ParamMap map(params, blk, 5, stretch_cols);

            ParamTablePtr itable = ParamTable::create();
            h->addAuxTable("stretch_fbhw_interaction", itable);
            itable->addProp("group1", IntType);
            itable->addProp("group2", IntType);
            itable->addProp("param",  IntType);

            const Json& g1 = blk.get("ffio_group1");
            const Json& g2 = blk.get("ffio_group2");
            
            int i,n = blk.get("__size__").as_int();
            for (i=0; i<n; i++) {
                Id itc = itable->addParam();
                Id p = map.add(i);
                itable->value(itc,"group1") = max_group;
                itable->value(itc,"group2") = max_group+1;
                itable->value(itc,"param" ) = p;
                IdList group1 = parse_ids(g1.elem(i).as_string());
                IdList group2 = parse_ids(g2.elem(i).as_string());
                BOOST_FOREACH(Id id, group1) {
                    Id t = table->addTerm(IdList(1,id), p);
                    table->termPropValue(t,"group") = max_group;
                }
                BOOST_FOREACH(Id id, group2) {
                    Id t = table->addTerm(IdList(1,id), p);
                    table->termPropValue(t,"group") = max_group+1;
                }
                max_group += 2;
            }
        }

    };

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
    RegisterFfio<Stretch>   _4("ffio_stretch_fbhw");
}

