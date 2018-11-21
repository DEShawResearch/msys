#include "../maeatoms.hxx"

namespace {

    struct Dihedral6 : public Ffio {

        void apply( SystemPtr h,
                    const Json& blk,
                    const SiteMap& sitemap,
                    const VdwMap& ) const {

            TermTablePtr table = AddTable(h,"dihedral6_trig");
            ParamMap map(table->params(), blk);
            MaeAtoms atoms(blk);
            
            const Json& fn = blk.get("ffio_funct");
            int i,n = blk.get("__size__").as_int();
            for (i=0; i<n; i++) {
                std::string f = fn.elem(i).as_string();
                to_lower(f);
                if (f!="trig") {
                    FFIO_ERROR("Expected ffio_funct='trig' in ffio_dihedrals6atom; got " << f);
                }
                sitemap.addUnrolledTerms( table, map.add(i), atoms.ids(i));
            }
        }
    };

    RegisterFfio<Dihedral6> _("ffio_dihedrals6atom");
}

