#include "../maeatoms.hxx"

namespace {

    //const char * maecols[] = { "ffio_a", "ffio_b", "ffio_cutoff" };
    struct Pseudopol : public Ffio {

        void apply( SystemPtr h,
                    const Json& blk,
                    const SiteMap& sitemap,
                    const VdwMap& ) const {

            TermTablePtr table = AddTable(h,"pseudopol_fermi");
            ParamMap map(table->params(), blk /* , 3, maecols */);
            MaeAtoms atoms(blk);
            
            const Json& fn = blk.get("ffio_funct");
            int i,n = blk.get("__size__").as_int();
            for (i=0; i<n; i++) {
                std::string f = fn.elem(i).as_string();
                boost::to_lower(f);
                if (f!="fermi") {
                    FFIO_ERROR("Expected ffio_funct='fermi' in ffio_pseudo_polarization; got " << f);
                }
                sitemap.addUnrolledTerms( table, map.add(i), atoms.ids(i));
            }
        }
    };

    RegisterFfio<Pseudopol> _("ffio_pseudo_polarization");
}

