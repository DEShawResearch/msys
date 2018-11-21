#ifndef desres_msys_hbond_hxx
#define desres_msys_hbond_hxx

#include "geom.hxx"

namespace desres { namespace msys { 

    struct HydrogenBond {
        double r;   /* donor-acceptor distance */
        double p;   /* angle formed by donor, hydrogen, and acceptor */

        /* angular deviations of H from bisector of lone-pair orbital */
        double to;  /* within-plane deviation */
        double ti;  /* out-of-plane deviation */

        HydrogenBond() : r(), p(), to(), ti() {}

        HydrogenBond(const double* dpos,  /* donor */
                 const double* apos,  /* acceptor */
                 const double* hpos,  /* hydrogen */
                 const double* cpos=NULL,  /* parent of acceptor */
                 const double* capos=NULL) /* grandparent of acceptor */
        : r(), p(), to(), ti()
        {
            if (!dpos || !apos) return;
            r = calc_distance(dpos, apos);

            if (hpos) {
                p = calc_angle(dpos, hpos, apos);
            }

            if (cpos && capos) {
                double proj[3];
                calc_projection(apos, cpos, capos, hpos, proj);
                ti = fabs(M_PI-calc_angle(proj, apos, cpos));
                to = calc_angle(hpos, apos, proj);
            }
        }

        double energy() const {

            /* distance term: Er = C/r^8 + D/r^6 */
            static const double Em =  2.8;
            static const double rm =  3.0;
            static const double rm2= rm*rm;
            static const double C =  3.0 * Em * pow(rm,8);
            static const double D =  4.0 * Em * pow(rm,6);

            double Er = -Em;
            double r2 = r * r;
            if (r2>=rm2) {
                double ir2 = 1.0/r2;
                double ir6 = ir2*ir2*ir2;
                double ir8 = ir6*ir2;
                Er = C*ir8 - D*ir6;
            }

            /* first directionality term */
            double Ep = p < M_PI/2 ? 0 : cos(p);
            Ep *= Ep;
            
            /* second directionality term */
            double Et = 0;
            if (ti<110 * M_PI/180) {
                double costo = cos(to);
                if (ti<90 * M_PI/180) {
                    Et = (0.9 + 0.1*sin(2*ti)) * costo;
                } else {
                    static const double K1 = 0.9/pow(cos(110*M_PI/180),6);
                    static const double K2 =     pow(cos(110*M_PI/180),2);

                    double costi = cos(ti);
                    Et = K2 - costi * costi;
                    Et *= Et*Et * K1 * costo;
                }
            }
            return Er * Ep * Et;
        }
    };
}}

#endif

