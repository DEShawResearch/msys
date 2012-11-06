#include "aromatic.hxx"
#include <stdexcept>
#include <list>
#include <cmath>

//#include <profiler/profiler.hxx>
#include "eigensystem.hxx"
#include "../elements.hxx"
#include "filtered_bonds.hxx"

namespace desres { namespace msys {

    AromaticAtom::Type
    AromaticAtom::Classify(int nb, int a0, int b0, int b1, int be) {
        if(nb<4){
            int vsum=a0-(3-nb);
            int bsum=b0+b1-2;
            int ebsum= be ? be-1 : be;
            if(vsum < 0 || vsum > 1 || bsum < 0 || bsum > 1 || ebsum>1 || (vsum==1 && bsum == 1)){
                // Not part of aromatic ring
                return INVALID;
            }
            if(ebsum==1 && vsum!=0 && bsum!=0) {
                std::stringstream msg;
                msg << "Bad combined vsum,bsum,ebsum in aromatic detection " << 
                    ": vsum=" << vsum << " bsum=" << bsum << " ebsum=" << ebsum << 
                    " b0=" << b0 << "  nb=" << nb << "  a0=" << a0 << " be="<< be << "  b1=" << b1;
                throw std::runtime_error(msg.str());
            }
            if(vsum==1){
                return X_TYPE;    // vsum=1, bsum=0, ebsum=0
            }else if(bsum==1){
                return Y_TYPE;    // vsum=0, bsum=1, ebsum=0
            }else if(ebsum==1){
                return YEXT_TYPE; // vsum=0, bsum=0, ebsum=1
            }else{
                return Z_TYPE;    // vsum=0, bsum=0, ebsum=0
            }
        }
        /* FIXME: This ends up excluding some some thiazole dioxide and isothazole dioxide compounds */
        return INVALID;
    }

    AromaticRing::Type
    AromaticRing::Classify(int nX, int nY, int nYe, int nZ) {
        /* cant be aromatic or antiaromatic without paired external electrons */
        if(nYe%2 == 1) {
            return NONAROMATIC;
        }

        /* Number of "extra" electrons in bonds around ring. MUST be even or something bad happened */ 
        if(nY%2 == 1) {
            std::stringstream msg;
            msg << "nY must be even in aromatic detection: nY = " << nY;
            throw std::runtime_error(msg.str());
        }

        int count= nX+(nY+nYe)/2;
        // Huckel 2n+1 rule for number of electron PAIRS, 'count-1' is even for aromatic
        if((count-1)%2 == 1){
            return ANTIAROMATIC;
        }
        return AROMATIC;       
    }

    bool ClassifyAromaticAtoms(SystemPtr mol, IdList const& atoms, 
                               int& nx, int& ny, int& nyext, int& nz) {

        nx=0;
        ny=0;
        nyext=0;
        nz=0;

        std::vector<int> typecounts(AromaticAtom::INVALID);
        IdList ringAtoms(atoms);
        size_t natoms= ringAtoms.size();
        Id previous= ringAtoms[natoms-1];
        /* simple detection for closed ring cycle in provided atoms list */
        if(ringAtoms[0]==ringAtoms[natoms-1]){
            natoms-=1;
            previous=ringAtoms[natoms-1];
        }else{
            ringAtoms.push_back(ringAtoms[0]);
        }
        for(size_t iatom=0; iatom<natoms; ++iatom){
            
            Id current = ringAtoms[iatom];
            Id next =  ringAtoms[iatom+1];
            
            IdList bonds=filteredBondsForAtom(mol,current);
            Id nb=bonds.size();

            atom_t & atm=mol->atom(current);
            int a0=DataForElement(atm.atomic_number).nValence-atm.formal_charge;
            unsigned b0,b1,be;
            b0=b1=be=0;
            BOOST_FOREACH( Id bid, bonds){
                bond_t & bond=mol->bond(bid);
                a0-=bond.order;
                Id other=bond.other(current);
                if(other==previous){
                    b0=bond.order;
                }else if(other==next){
                    b1=bond.order;
                }else if(nb==3 && mol->atom(current).atomic_number==6 &&
                         mol->atom(other).atomic_number == 6){
                    be=bond.order;
                }
            }
            assert(b0!=0 && b1!=0);
            assert(a0>=0 && a0%2==0);
            a0/=2;

            AromaticAtom::Type atype=AromaticAtom::Classify(nb,a0,b0,b1,be);
            if(atype>=AromaticAtom::INVALID){
                return false;
            }
            typecounts[atype]++;

            previous=current;
        }

        nx = typecounts[AromaticAtom::X_TYPE];
        ny = typecounts[AromaticAtom::Y_TYPE];
        nyext = typecounts[AromaticAtom::YEXT_TYPE];
        nz = typecounts[AromaticAtom::Z_TYPE];
        return true;
    }

    double ComputeRingPlanarity(const SystemPtr mol, const IdList& aids){
        //static desres::profiler::Symbol _("::ring_planarity_descriptor");
        //desres::profiler::Clock __(_);

        size_t nids=aids.size();

        /* simple detection for closed ring cycle in provided atoms list */
        if(nids>1 && aids[0]==aids[nids-1]) nids--;
        if(nids<3) return 0.0;
 
        double xctr=0.0,yctr=0.0,zctr=0.0;
        for(size_t idx=0; idx<nids;++idx){
            const atom_t& atm=mol->atom(aids[idx]);
            xctr+=atm.x; yctr+=atm.y; zctr+=atm.z;
        }
        xctr/=nids; yctr/=nids; zctr/=nids;
        double I[9]={0,0,0,0,0,0,0,0,0};
        for(size_t idx=0; idx<nids;++idx){
            const atom_t& atm=mol->atom(aids[idx]);
            double x=atm.x-xctr;
            double y=atm.y-yctr;
            double z=atm.z-zctr;
            I[0] += (y*y + z*z);
            I[1] += -(x*y);
            I[2] += -(x*z);
            I[4] += (x*x + z*z);
            I[5] += -(y*z);
            I[8] += (x*x + y*y);
        }
        I[3]=I[1];
        I[6]=I[2];
        I[7]=I[5];
        
        double v[3];
        real_symmetric_eigenvalues_3x3(I,v,NULL,NULL);
        // perpendicular axis theorem
        double planar=fabs(v[0]-(v[1]+v[2]));
        //printf("planarity descriptor for %zu atoms = %f  (Eigenvalues= %f %f %f for %zu atoms)\n",nids,planar,v[0],v[1],v[2],nids);
        return planar;
    }

    AromaticRing::Type ClassifyAromaticRing(SystemPtr mol, IdList const& atms) {
        int nx, ny, nyext, nz;
        if (ClassifyAromaticAtoms(mol, atms, nx, ny, nyext, nz)) {
            return AromaticRing::Classify(nx, ny, nyext, nz);
        }
        return AromaticRing::NONAROMATIC;
    }
}}
