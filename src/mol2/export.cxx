#include "../mol2.hxx"
#include "elements.hxx"
#include <boost/foreach.hpp>
#include <stdio.h>
#include <math.h>

using namespace desres::msys;

static 
double calc_dot_prod( const double* A, const double* B ) {
    return A[0]*B[0] + A[1]*B[1] + A[2]*B[2];
}

static
void calc_cross_prod( double* C, const double* A, const double* B ) {
    double cx, cy, cz;
    cx =  A[1]*B[2] - B[1]*A[2];
    cy = -A[0]*B[2] + B[0]*A[2];
    cz =  A[0]*B[1] - B[0]*A[1];
    C[0] = cx;
    C[1] = cy;
    C[2] = cz;
}

static
double calc_vec_angle( const double* r1, const double* r2 ) {
    double psin, pcos;

    double r3[3];
    calc_cross_prod( r3, r1, r2  );
    psin = sqrt(calc_dot_prod(r3,r3));
    pcos = calc_dot_prod( r1, r2 );
    return atan2(psin,pcos);
}

double calc_angle( const double* A, const double* B, const double* C ) {
    int i;
    double r1[3], r2[3];
    for (i=0; i<3; i++) {
        r1[i]=A[i]-B[i];
        r2[i]=C[i]-B[i];
    }
    return calc_vec_angle( r1, r2 );
}


/* reference: 
 * http://www.sdsc.edu/CCMS/Packages/cambridge/pluto/atom_types.html
 */
const char* guess_atom_type(SystemPtr mol, Id id) {
    atom_t const& atm = mol->atom(id);
    residue_t const& res = mol->residue(atm.residue);
    Id nbnd=mol->bondCountForAtom(id);
    Id nres=0; /* number of bonds with resonant_order 1.5 */
    Id nsng=0;  /* number of single bonds */
    Id ndbl=0;  /* number of double bonds */
    Id ntrp=0;  /* number of triple bonds */
    Id nnit=0;  /* number of bonds to nitrogen */
    Id ncar=0;  /* number of bonds to carbon */
    Id npho=0;  /* number of bonds to phosphorus */
    Id nhyd=0;  /* number of bonds to hydrogen */
    BOOST_FOREACH(Id bid, mol->bondsForAtom(id)) {
        bond_t const& bond = mol->bond(bid);
        if (bond.resonant_order==1.5) ++nres;
        if (bond.order==1) ++nsng;
        if (bond.order==2) ++ndbl;
        if (bond.order==3) ++ntrp;
        atom_t const& other = mol->atom(bond.other(id));
        if      (other.atomic_number==1) ++nhyd;
        else if (other.atomic_number==6) ++ncar;
        else if (other.atomic_number==7) ++nnit;
        else if (other.atomic_number==15) ++npho;
    }
    //if (atm.atomic_number>1) {
        //printf("atom %4u %8s bonds=%u res=%u sng=%u dbl=%u\n",
                //id, atm.name.c_str(), nbnd, nres, nsng, ndbl);
    //}
    switch (atm.atomic_number) {
        case 0:
            return "Du";    /* 1.1 */
        case 1:
            return "H";     /* 1.3 */
        case 15:
            return "P.3";   /* 1.4 */
        case 27:
            return "Co.oh"; /* 1.5 */
        case 44:
            return "Ru.oh";
        case 6:             /* 1.6 */
            if (nbnd>=4 && nsng==nbnd) return "C.3";    /* 1.6.1 */
            if (nbnd==3 && nres==0 && nnit==3) {        /* 1.6.2 */
                /* check that each nitrogen forms bonds to two other atoms,
                 * neither of which is oxygen */
                bool ok=true;
                BOOST_FOREACH(Id atm, mol->bondedAtoms(id)) {
                    if (mol->bondCountForAtom(atm)!=3) { ok=false; break; }
                    BOOST_FOREACH(Id nbr, mol->bondedAtoms(atm)) {
                        if (mol->atom(nbr).atomic_number==8) {ok=false; break;}
                    }
                }
                if (ok) return "C.cat";
            }
            if (nbnd>=2 && nres==2) return "C.ar";
            if ((nbnd==1 || nbnd==2) && ntrp==1) return "C.1";
            return "C.2";
        case 8:             /* 1.7 */
            if (res.type==ResidueWater) return "O.t3p";
            if (nbnd==1) {  /* 1.7.1 */
                if (ncar==1 || npho==1) {       /* 1.7.1.1, 1.7.1.2 */
                    Id parent = mol->bondedAtoms(id).at(0);
                    if (mol->bondCountForAtom(parent)==3) {
                        bool found=false;
                        BOOST_FOREACH(Id nbr, mol->bondedAtoms(parent)) {
                            if (nbr==id) continue;
                            if (mol->atom(nbr).atomic_number==8 &&
                                mol->bondCountForAtom(nbr)==1) {
                                found = !found; /* want to find exactly 1 */
                            }
                        }
                        if (found) return "O.co2";
                    }
                }
            }
            if (nbnd>=2 && nsng==nbnd) return "O.3";    /* 1.7.2 */
            return "O.2";                               /* 1.7.3 */
        case 7:
            if (nbnd==4 && nsng==nbnd) return "N.4";    /* 1.8.1. */
            if (nres==2) return "N.ar";                 /* 1.8.2 */
            if (nbnd==1 && ntrp==1) return "N.1";       /* 1.8.3 */
            if (nbnd==2 && (ndbl==2 || (nsng==1 && ntrp==1))) {
                return "N.1";                           /* 1.8.4 */
            }
            if (nbnd==3 && nhyd<=1 && ncar>=1) {        /* 1.8.5 */
                /* check for amide */
                int namide=0;
                BOOST_FOREACH(Id car, mol->bondedAtoms(id)) {
                    if (mol->atom(car).atomic_number==6) {
                        BOOST_FOREACH(Id bnd, mol->bondsForAtom(car)) {
                            Id other = mol->bond(bnd).other(car);
                            if (mol->bond(bnd).order==2 && (
                                mol->atom(other).atomic_number==8 ||
                                mol->atom(other).atomic_number==16)) {
                                ++namide;
                            }
                        }
                    }
                }
                if (namide==1) return "N.am";
            }
            if (nbnd==3) {                              /* 1.8.6 */
                if (nsng==2) return "N.pl3";            /* 1.8.6.1 */
                if (nsng==3) {                          /* 1.8.6.2 */
                    int nh=0, nd=0;
                    BOOST_FOREACH(Id nbr, mol->bondedAtoms(id)) {
                        if (mol->atom(nbr).atomic_number==1) ++nh;
                        else BOOST_FOREACH(Id bnd, mol->bondsForAtom(nbr)) {
                            if (mol->bond(bnd).resonant_order>1) {
                                ++nd;
                                break;
                            }
                        }
                    }
                    if (nh>=1 && nd==1) return "N.pl3"; /* 1.8.6.2.1 */
                    if (nh==0 && nd==1) {               /* 1.8.6.2.2 */
                        double tot=0;
                        IdList nbrs = mol->bondedAtoms(id);
                        const double* A = &mol->atom(nbrs[0]).x;
                        const double* B = &mol->atom(nbrs[1]).x;
                        const double* C = &mol->atom(nbrs[2]).x;
                        const double* o = &mol->atom(id).x;
                        tot += calc_angle(A,o,B);
                        tot += calc_angle(B,o,C);
                        tot += calc_angle(C,o,A);
                        tot *= 180/M_PI;
                        if (tot>=350) return "N.pl3";
                    }
                }
                return "N.3";                           /* 1.8.6.3 */
            }
            return "N.2";                               /* 1.8.7 */
        case 16:                                        /* 1.9 */
            {
                int nox=0;
                BOOST_FOREACH(Id nbr, mol->bondedAtoms(id)) {
                    if (mol->atom(nbr).atomic_number==8 &&
                        mol->bondCountForAtom(nbr)==1) ++nox;
                }
                if (nbnd==3 && nox==1) return "S.o";    /* 1.9.1 */
                if (nbnd==4 && nox==2) return "S.o2";   /* 1.9.2 */
                if (nbnd>=2 && nsng==nbnd) return "S.3"; /* 1.9.3 */
                return "S.2";                           /* 1.9.4 */
            }
        default: ;
    }
    return AbbreviationForElement(atm.atomic_number);
}

static const char* guess_bond_type( SystemPtr mol, Id bnd,
                                    std::string const& itype, 
                                    std::string const& jtype) {

    if (itype=="N.am" || jtype=="N.am") {
        Id n = mol->bond(bnd).i;
        if (mol->atom(n).atomic_number!=7) n = mol->bond(bnd).j;
        Id car = mol->bond(bnd).other(n);
        if (mol->atom(car).atomic_number==6) {
            BOOST_FOREACH(Id bnd, mol->bondsForAtom(car)) {
                Id other = mol->bond(bnd).other(car);
                if (mol->bond(bnd).order==2 && (
                    mol->atom(other).atomic_number==8 ||
                    mol->atom(other).atomic_number==16)) {
                    return "am";
                }
            }
        }
    }
    if (itype=="C.ar" && jtype=="C.ar") return "ar";
    int order = mol->bond(bnd).order;
    if (order==1) return "1";
    if (order==2) return "2";
    if (order==3) return "3";
    return "un";
}

void desres::msys::ExportMol2( SystemPtr mol, std::string const& path,
                               Provenance const& provenance) {

    FILE* fd = fopen(path.c_str(), "wb");
    if (!fd) MSYS_FAIL("Could not open '" << "' for writing.");
    boost::shared_ptr<FILE> dtor(fd, fclose);

    /* mapping of atom id to atom type */
    std::vector<const char*> atypes(mol->maxAtomId());

    /* molecule record */
    fprintf(fd, "@<TRIPOS>MOLECULE\n");
    fprintf(fd, "%s\n", mol->name.c_str());
    fprintf(fd, " %4d %4d %4d\n", 
            mol->atomCount(), mol->bondCount(), mol->residueCount());
    fprintf(fd, "%s\n", mol->residueCount()==1 ? "SMALL" : "BIOPOLYMER");
    fprintf(fd, "USER_CHARGES\n");
    fprintf(fd, "%s\n", "****");    /* status bits */
    fprintf(fd, "\n");              /* comment */

    /* atom records */
    IdList idmap(mol->maxAtomId(), 0);
    Id index=0;

    fprintf(fd, "@<TRIPOS>ATOM\n");
    for (Id i=0; i<mol->maxAtomId(); i++) {
        if (!mol->hasAtom(i)) continue;
        idmap[i] = ++index;
        const atom_t& atm = mol->atom(i);
        const residue_t& res = mol->residue(atm.residue);

        /* Use atom name unless it's invalid */
        std::string aname(atm.name);
        if (aname.size()<1 || aname.size()>7 || !isalpha(aname[0])) { 
            std::stringstream ss;
            ss << AbbreviationForElement(atm.atomic_number) << i;
            aname = ss.str();
        }

        /* Use residue name unless it's invalid */
        std::string rname(res.name);
        if (rname.size()<1 || rname.size()>7 || !isalpha(rname[0])) {
            rname = "UNK";
        }

        /* guess an atom type */
        const char* type = guess_atom_type(mol, i);
        atypes.at(i)=type;

        /* write the atom line */
        fprintf(fd, 
           "%7d %-4s      %8.4f  %8.4f  %8.4f %-6s %4d  %4s  %8.4f\n", 
           index, aname.c_str(), atm.x, atm.y, atm.z,
           type, atm.residue+1, rname.c_str(), atm.charge);
    }

    /* bond records */
    fprintf(fd, "@<TRIPOS>BOND\n");
    for (Id i=0; i<mol->maxBondId(); i++) {
        if (!mol->hasBond(i)) continue;
        bond_t const& bnd = mol->bond(i);
        Id ai = idmap[bnd.i];
        Id aj = idmap[bnd.j];
        const char* itype = atypes.at(bnd.i);
        const char* jtype = atypes.at(bnd.j);
        const char* btype = guess_bond_type(mol, i, itype, jtype);

        fprintf(fd, "%5u %5u %5u %s\n", i+1, ai, aj, btype);
    }

    /* substructure */
    fprintf(fd, "@<TRIPOS>SUBSTRUCTURE\n");
    for (Id i=0; i<mol->maxResidueId(); i++) {
        if (!mol->atomCountForResidue(i)) continue;
        residue_t const& res = mol->residue(i);
        fprintf(fd, "%7u %-4s %7u %-8s\n",
                i+1,                                /* residue id */
                res.name.c_str(),                   /* residue name */
                mol->atomsForResidue(i).at(0)+1,    /* root atom */
                mol->residueCount()==1 ? "GROUP" : "RESIDUE");
    }
}

