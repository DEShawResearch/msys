#include "../maeatoms.hxx"

namespace {

    struct Constraints : public Ffio {

        void apply( SystemPtr h,
                    const Json& blk,
                    const SiteMap& sitemap,
                    const VdwMap& ) const {

            MaeAtoms atoms(blk);

            typedef boost::shared_ptr<ParamMap> ParamMapPtr;
            typedef std::pair<TermTablePtr, ParamMapPtr> type_t;
            typedef std::map<std::string, type_t> TypeMap;
            TypeMap typemap;

            const Json& fn = blk.get("ffio_funct");
            int i,n = blk.get("__size__").as_int();
            for (i=0; i<n; i++) {
                std::string f = fn.elem(i).as_string();
                boost::to_lower(f);
                type_t type = typemap[f];
                if (!type.first) {
                    std::string name = std::string("virtual_")+f;
                    type.first = AddTable(h, name);
                    ParamTablePtr params = type.first->params();
                    type.second = ParamMapPtr(new ParamMap(
                                params, blk, params->propCount()));
                    typemap[f] = type;
                }
                Id m = type.first->atomCount();
                Id A = type.second->add(i);
                const IdList& ids = atoms.ids(i,m);
                sitemap.addUnrolledTerms( type.first, A, ids);
                /* put the pseudo second since it has higher id */
                sitemap.addUnrolledPseudoBonds( h, ids[1], ids[0] );
            }
        }
    };

    RegisterFfio<Constraints> _("ffio_virtuals");
}

