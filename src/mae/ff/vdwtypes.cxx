#include "../ff.hxx"
#include <cstdio>

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
                funct="vdw_exp6";
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

#if 0
            if (vdwmap.combined().size()) {
                /* create nonbonded_combined_param extra table with columns
                 * param1, param2, properties... */
                ent::DictPtr d;
                static const char table_name[] = "nonbonded_combined_param";
                if (h.hasExtraTable(table_name)) {
                    d=h.getExtraTable(table_name);
                } else {
                    d.reset(new ent::Dict);
                    h.addExtraTable(table_name, d);
                    d->addColumn( "param1", ent::DICT_VALUE_INT );
                    d->addColumn( "param2", ent::DICT_VALUE_INT );

                    for (int i=0; i<info->nprops; i++) {
                        d->addColumn(info->propnames[i], ent::DICT_VALUE_FLOAT);
                    }
                }
                EntryMap::const_iterator i, j, e=emap.end();
                for (i=emap.begin(); i!=e; ++i) {
                    const VdwType& itype = i->first;
                    int iid = nb.getEntryId(i->second);
                    for (j=emap.begin(); j!=e; ++j) {
                        const VdwType& jtype = j->first;
                        int jid = nb.getEntryId(j->second);
                        if (vdwmap.has_combined(itype,jtype)) {
                            const VdwParam& p = vdwmap.param(itype,jtype);
                            ent::DictEntry& row = *(d->addRow());
                            row[0] = ent::DictValue(iid);
                            row[1] = ent::DictValue(jid);
                            for (int iprop=0; iprop<info->nprops; iprop++) {
                                row[2+iprop]=ent::DictValue(p.at(iprop));
                            }
                        }
                    }
                }
            }
#endif
        }
    };
    RegisterFfio<VdwTypes> _("ffio_vdwtypes");
    RegisterFfio<VdwTypesCombined> _2("ffio_vdwtypes_combined");
}

