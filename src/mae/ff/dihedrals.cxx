#include "../ff.hxx"

namespace {

    void convert_opls( std::vector<double>& v ) {
        assert(v.size()==5);
        double phi = v[0];
        double c1  = v[1];
        double c2  = v[2];
        double c3  = v[3];
        double c4  = v[4];

        v.resize(8);

        v[0] = phi;
        v[1] =  0.5*( c1+c2+c3+c4 );
        v[2] =  0.5*c1;
        v[3] = -0.5*c2;
        v[4] =  0.5*c3;
        v[5] = -0.5*c4;
        v[6] =  0.0;
        v[7] =  0.0;
    }

    struct Dihedrals : public Ffio {

        void apply( SystemPtr h,
                    const Json& blk,
                    const SiteMap& sitemap,
                    const VdwMap&, bool alchemical  ) const {

            static const char * maecols[] = {
                "ffio_c0", "ffio_c1", "ffio_c2", "ffio_c3",
                "ffio_c4", "ffio_c5", "ffio_c6", "ffio_c7"
            };

            TermTablePtr dihedrals, impropers, anharms, opls;
            ParamMapPtr dmap, imap, amap, omap;

            const Json& ai = blk.get("ffio_ai");
            const Json& aj = blk.get("ffio_aj");
            const Json& ak = blk.get("ffio_ak");
            const Json& al = blk.get("ffio_al");
            const Json& fn = blk.get("ffio_funct");
            IdList ids(4);

            TermTablePtr table;
            int i,n = blk.get("__size__").as_int();
            for (i=0; i<n; i++) {
                std::string f = fn.elem(i).as_string();
                boost::to_lower(f);
                Id A=BadId, B=BadId;

                if (f=="proper_trig" || f=="improper_trig") {
                    if (!dihedrals) {
                        dihedrals = AddTable(h,"dihedral_trig");
                        dmap.reset( 
                           new ParamMap(dihedrals->params(), blk, 8, maecols));
                    }
                    table=dihedrals;
                    A=dmap->add(i);
                    if (alchemical) B = dmap->addAlchemical(i);

                } else if (f=="improper_harm") {
                    if (!impropers) {
                        impropers = AddTable(h,"improper_harm");
                        imap.reset(
                            new ParamMap(impropers->params(), blk, 2, maecols));
                    }
                    table=impropers;
                    A=imap->add(i);
                    if (alchemical) B = imap->addAlchemical(i);

                } else if (f=="improper_anharm") {
                    if (!anharms) {
                        anharms = AddTable(h,"improper_anharm");
                        amap.reset(
                            new ParamMap(anharms->params(),   blk, 2, maecols));
                    }
                    table=anharms;
                    A=amap->add(i);
                    if (alchemical) B = amap->addAlchemical(i);

                } else if (f=="opls_proper" || f=="opls_improper") {
                    if (!opls) {
                        opls = AddTable(h,"dihedral_trig");
                        omap.reset(
                            new ParamMap(opls->params(),      blk, 5, maecols));
                    }
                    table=opls;
                    A=omap->add(i, convert_opls);
                    if (alchemical) B = omap->addAlchemical(i, convert_opls);
                } else {
                    FFIO_ERROR( "Unsupported ffio_funct '" << f 
                                << "' in ffio_dihedrals");
                }
                ids[0]=ai.elem(i).as_int();
                ids[1]=aj.elem(i).as_int();
                ids[2]=ak.elem(i).as_int();
                ids[3]=al.elem(i).as_int();
                sitemap.addUnrolledTerms( table, A, ids, false, B );
            }
        }
    };

    RegisterFfio<Dihedrals> _("ffio_dihedrals");
}

