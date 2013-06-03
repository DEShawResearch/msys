#include "../ff.hxx"
#include "../../override.hxx"
//#include <cstdio>
#include <boost/foreach.hpp>
#include <assert.h>

namespace {

    struct VdwTypesCombined : public Ffio {
        void apply( SystemPtr h,
                    const Json& blk,
                    const SiteMap& sitemap,
                    const VdwMap&  vdwmap ) const {}
    };

    struct VdwTypes : public Ffio {
        void apply( SystemPtr h, 
                    const Json& blk,
                    const SiteMap& sitemap,
                    const VdwMap&  vdwmap ) const {

            std::string funct;
            if (vdwmap.funct()=="lj12_6_sig_epsilon") {
                funct="vdw_12_6";
            } else if (vdwmap.funct()=="exp_6x") {
                funct="vdw_exp_6";
            } else if (vdwmap.funct()=="polynomial_cij") {
                funct="polynomial_cij";
            } else {
                std::stringstream ss;
                ss << "Unrecognized mae vdw_funct '" << vdwmap.funct()
                    << "'";
                throw std::runtime_error(ss.str());
            }

            TermTablePtr table = AddNonbonded(h,funct, vdwmap.rule());
            TermTablePtr atable;
            ParamTablePtr params = table->params();
            Id typecol = params->addProp("type", StringType);

            typedef std::map<VdwType, Id> TypeMap;
            TypeMap map;

            for (Id i=0; i<sitemap.nsites(); i++) {
                Id site = i+1;
                const VdwType* types[2];
                types[0] = &vdwmap.type(site);
                types[1] = &vdwmap.typeB(site);
                Id p[2]={BadId,BadId};
                for (int j=0; j<2; j++) {
                    const VdwType& type = *types[j];
                    if (type=="") continue;
                    TypeMap::const_iterator e = map.find(type);
                    if (e==map.end()) {
                        p[j] = map[type] = params->addParam();
                        params->value(p[j], typecol) = type;
                        const VdwParam& vals = vdwmap.param(type);
                        for (Id k=0; k<vals.size(); k++) {
                            params->value(p[j], k) = vals[k];
                        }
                    } else {
                        p[j] = e->second;
                    }
                }
                IdList ids(1,site);
                sitemap.addUnrolledTerms(table, p[0], ids );

                if (!bad(p[1])) {
                    /* alchemical nonbonded term */
                    if (!atable) {
                        atable = h->addTable(
                                "alchemical_nonbonded", 1, params);
                        atable->category = NONBONDED;
                        atable->addTermProp("chargeB", FloatType);
                        atable->addTermProp("moiety", IntType);
                    }
                    IdList realids(1,sitemap.site(site));
                    Id term = atable->addTerm(realids, p[1]);
                    atable->termPropValue(term,"chargeB") = 
                        vdwmap.chargeB(site);

                }
            }

            if (vdwmap.combined().size()) {
                OverrideTablePtr o = table->overrides();
                ParamTablePtr p = o->params();
                for (Id i=0; i<params->propCount(); i++) {
                    p->addProp(params->propName(i), params->propType(i));
                }
                BOOST_FOREACH(TypeMap::value_type ti, map) {
                    VdwType const& itype = ti.first;
                    BOOST_FOREACH(TypeMap::value_type tj, map) {
                        VdwType const& jtype = tj.first;
                        if (vdwmap.has_combined(itype, jtype)) {
                            Id row = p->addParam();
                            VdwParam const& props = vdwmap.param(itype, jtype);
                            for (Id i=0; i<props.size(); i++) {
                                p->value(row,i)=props[i];
                            }
                            Id pi = ti.second;
                            Id pj = tj.second;
                            o->set(IdPair(pi,pj), row);
                        }
                    }
                }
            }
        }
    };
    RegisterFfio<VdwTypes> _("ffio_vdwtypes");
    RegisterFfio<VdwTypesCombined> _2("ffio_vdwtypes_combined");
}

