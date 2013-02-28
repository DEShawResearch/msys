#include "../ff.hxx"

namespace {

    void convert_opls( double* v ) {
        double phi = v[0];
        double c1  = v[1];
        double c2  = v[2];
        double c3  = v[3];
        double c4  = v[4];

        v[0] = phi;
        v[1] =  0.5*( c1+c2+c3+c4 );
        v[2] =  0.5*c1;
        v[3] = -0.5*c2;
        v[4] =  0.5*c3;
        v[5] = -0.5*c4;
        v[6] =  0.0;
        v[7] =  0.0;
    }

    template <bool alchemical>
    struct Dihedrals : public Ffio {

        typedef std::vector<const char *> StrList;
        static StrList make_columns(int n) {
            static const char * regcols[] = {
                "ffio_c0", "ffio_c1", "ffio_c2", "ffio_c3",
                "ffio_c4", "ffio_c5", "ffio_c6", "ffio_c7"
            };
            static const char * alccols[] = {
                "ffio_c0B", "ffio_c1B", "ffio_c2B", "ffio_c3B",
                "ffio_c4B", "ffio_c5B", "ffio_c6B", "ffio_c7B"
            };
            StrList cols(regcols,regcols+n);
            if (alchemical) cols.insert(cols.end(), alccols, alccols+n);
            return cols;
        }

        void apply( SystemPtr h,
                    const Json& blk,
                    const SiteMap& sitemap,
                    const VdwMap& ) const {


            TermTablePtr dihedrals, impropers, anharms;
            ParamMapPtr dmap, imap, amap;

            const Json& ai = blk.get("ffio_ai");
            const Json& aj = blk.get("ffio_aj");
            const Json& ak = blk.get("ffio_ak");
            const Json& al = blk.get("ffio_al");
            const Json& fn = blk.get("ffio_funct");
            IdList ids(4);

            std::string prefix;
            if (alchemical) prefix = "alchemical_";

            StrList oplscols = make_columns(5);

            TermTablePtr table;
            int i,n = blk.get("__size__").as_int();
            for (i=0; i<n; i++) {
                std::string f = fn.elem(i).as_string();
                boost::to_lower(f);
                Id A=BadId;

                if (f=="proper_trig" || f=="improper_trig") {
                    if (!dihedrals) {
                        dihedrals = AddTable(h,prefix+"dihedral_trig");
                        StrList cols = make_columns(8);
                        dmap.reset( 
                           new ParamMap(dihedrals->params(), blk, 
                               cols.size(), &cols[0]));
                    }
                    table=dihedrals;
                    A=dmap->add(i);

                } else if (f=="proper_harm" || f=="improper_harm") {
                    if (!impropers) {
                        impropers = AddTable(h,prefix+"improper_harm");
                        StrList cols = make_columns(2);
                        imap.reset(
                            new ParamMap(impropers->params(), blk, 
                                cols.size(), &cols[0]));
                    }
                    table=impropers;
                    A=imap->add(i);

                } else if (f=="improper_anharm") {
                    if (!anharms) {
                        anharms = AddTable(h,prefix+"improper_anharm");
                        StrList cols = make_columns(2);
                        amap.reset(
                            new ParamMap(anharms->params(), blk,
                                cols.size(), &cols[0]));
                    }
                    table=anharms;
                    A=amap->add(i);

                } else if (f=="opls_proper" || f=="opls_improper") {
                    if (!dihedrals) {
                        dihedrals = AddTable(h,prefix+"dihedral_trig");
                        StrList cols = make_columns(8);
                        dmap.reset(
                            new ParamMap(dihedrals->params(), blk,
                                cols.size(), &cols[0]));
                    }
                    table=dihedrals;
                    std::vector<double> p(alchemical ? 16 : 8);
                    for (int j=0; j<5; j++) {
                        p[j] = blk.get(oplscols[j]).elem(i).as_float();
                        if (alchemical) {
                            p[j+8] = blk.get(oplscols[5+j]).elem(i).as_float();
                        }
                    }
                    convert_opls(&p[0]);
                    if (alchemical) convert_opls(&p[8]);
                    A=dmap->add(p);
                } else {
                    FFIO_ERROR( "Unsupported ffio_funct '" << f 
                                << "' in ffio_dihedrals");
                }
                ids[0]=ai.elem(i).as_int();
                ids[1]=aj.elem(i).as_int();
                ids[2]=ak.elem(i).as_int();
                ids[3]=al.elem(i).as_int();
                sitemap.addUnrolledTerms( table, A, ids );
            }
        }
    };

    RegisterFfio<Dihedrals<false> > _1("ffio_dihedrals");
    RegisterFfio<Dihedrals<true> > _2("ffio_dihedrals_alchemical");
}

