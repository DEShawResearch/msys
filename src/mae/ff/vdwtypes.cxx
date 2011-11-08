#include "../ff.hxx"
//#include <cstdio>
#include <boost/foreach.hpp>

namespace {

    struct VdwTypesCombined : public Ffio {
        void apply( SystemPtr h,
                    const Json& blk,
                    const SiteMap& sitemap,
                    const VdwMap&  vdwmap, bool alchemical ) const {}
    };

    struct VdwTypes : public Ffio {
        void apply( SystemPtr h, 
                    const Json& blk,
                    const SiteMap& sitemap,
                    const VdwMap&  vdwmap, bool alchemical ) const {

            std::string funct;
            if (vdwmap.funct()=="lj12_6_sig_epsilon") {
                funct="vdw_12_6";
            } else if (vdwmap.funct()=="exp_6x") {
                funct="vdw_exp_6";
            } else if (vdwmap.funct()=="polynomial_pij") {
                funct="polynomial_cij";
            } else {
                std::stringstream ss;
                ss << "Unrecognized mae vdw_funct '" << vdwmap.funct()
                    << "'";
                throw std::runtime_error(ss.str());
            }

            TermTablePtr table = AddNonbonded(h,funct, vdwmap.rule());
            ParamTablePtr params = table->params();

            typedef std::map<VdwType, Id> TypeMap;
            TypeMap map;

            for (Id i=0; i<sitemap.nsites(); i++) {
                Id site = i+1;
                const VdwType* types[2];
                types[0] = &vdwmap.type(site);
                types[1] = &vdwmap.typeB(site);
                Id p[2];
                for (int j=0; j<2; j++) {
                    const VdwType& type = *types[j];
                    TypeMap::const_iterator e = map.find(type);
                    if (e==map.end()) {
                        p[j] = map[type] = params->addParam();
                        const VdwParam& vals = vdwmap.param(type);
                        for (Id k=0; k<params->propCount(); k++) {
                            params->value(p[j], k) = vals.at(k);
                        }
                    } else {
                        p[j] = e->second;
                    }
                }
                IdList ids(1,site);
                sitemap.addUnrolledTerms(table, p[0], ids, false, p[1] );
            }

            if (vdwmap.combined().size()) {
                /* create nonbonded_combined_param extra table with columns
                 * param1, param2, properties... */
                ParamTablePtr p = ParamTable::create();
                h->addExtra("nonbonded_combined_param", p);
                p->addProp("param1", IntType);
                p->addProp("param2", IntType);
                for (Id i=0; i<params->propCount(); i++) {
                    p->addProp(params->propName(i), params->propType(i));
                }
                BOOST_FOREACH(TypeMap::value_type ti, map) {
                    VdwType const& itype = ti.first;
                    BOOST_FOREACH(TypeMap::value_type tj, map) {
                        VdwType const& jtype = tj.first;
                        if (vdwmap.has_combined(itype, jtype)) {
                            Id row = p->addParam();
                            p->value(row,0)=ti.second;
                            p->value(row,1)=tj.second;
                            VdwParam const& props = vdwmap.param(itype, jtype);
                            for (Id i=0; i<props.size(); i++) {
                                p->value(row,2+i)=props[i];
                            }
                        }
                    }
                }
            }
        }
    };
    RegisterFfio<VdwTypes> _("ffio_vdwtypes");
    RegisterFfio<VdwTypesCombined> _2("ffio_vdwtypes_combined");
}

