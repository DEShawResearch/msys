
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>
#include "charmm_file.h"
#include "defs.hxx"

#include <boost/shared_ptr.hpp>

#include <stdexcept>

#if defined(_MSC_VER)
#define strcasecmp  stricmp
#define strncasecmp strnicmp
#endif

#define TOKLEN 100
#define BUFLEN 200

using namespace desres::msys::builder;

namespace {
  struct element {
    double daltons;
    const char* abbreviation;
    const char* name;
  };
}

// url = "http://physics.nist.gov/cgi-bin/Elements/elInfo.pl?element=%d&context=noframes"%element
static struct element amu[] = {
  {1.00794,"H","Hydrogen"},
  {4.002602,"He","Helium"},
  {6.941,"Li","Lithium"},
  {9.012182,"Be","Beryllium"},
  {10.811,"B","Boron"},
  {12.0107,"C","Carbon"},
  {14.0067,"N","Nitrogen"},
  {15.9994,"O","Oxygen"},
  {18.9984032,"F","Fluorine"},
  {20.1797,"Ne","Neon"},
  {22.989770,"Na","Sodium"},
  {24.3050,"Mg","Magnesium"},
  {26.981538,"Al","Aluminum"},
  {28.0855,"Si","Silicon"},
  {30.973761,"P","Phosphorus"},
  {32.065,"S","Sulfur"},
  {35.453,"Cl","Chlorine"},
  {39.0983,"K","Potassium"},
  {39.948,"Ar","Argon"},
  {40.078,"Ca","Calcium"},
  {44.955910,"Sc","Scandium"},
  {47.867,"Ti","Titanium"},
  {50.9415,"V","Vanadium"},
  {51.9961,"Cr","Chromium"},
  {54.938049,"Mn","Manganese"},
  {55.845,"Fe","Iron"},
  {58.6934,"Ni","Nickel"},
  {58.933200,"Co","Cobalt"},
  {63.546,"Cu","Copper"},
  {65.409,"Zn","Zinc"},
  {69.723,"Ga","Gallium"},
  {72.64,"Ge","Germanium"},
  {74.92160,"As","Arsenic"},
  {78.96,"Se","Selenium"},
  {79.904,"Br","Bromine"},
  {83.798,"Kr","Krypton"},
  {85.4678,"Rb","Rubidium"},
  {87.62,"Sr","Strontium"},
  {88.90585,"Y","Yttrium"},
  {91.224,"Zr","Zirconium"},
  {92.90638,"Nb","Niobium"},
  {95.94,"Mo","Molybdenum"},
  {101.07,"Ru","Ruthenium"},
  {102.90550,"Rh","Rhodium"},
  {106.42,"Pd","Palladium"},
  {107.8682,"Ag","Silver"},
  {112.411,"Cd","Cadmium"},
  {114.818,"In","Indium"},
  {118.710,"Sn","Tin"},
  {121.760,"Sb","Antimony"},
  {126.90447,"I","Iodine"},
  {127.60,"Te","Tellurium"},
  {131.293,"Xe","Xenon"},
  {132.90545,"Cs","Cesium"},
  {137.327,"Ba","Barium"},
  {138.9055,"La","Lanthanum"},
  {140.116,"Ce","Cerium"},
  {140.90765,"Pr","Praseodymium"},
  {144.24,"Nd","Neodymium"},
  {150.36,"Sm","Samarium"},
  {151.964,"Eu","Europium"},
  {157.25,"Gd","Gadolinium"},
  {158.92534,"Tb","Terbium"},
  {162.500,"Dy","Dysprosium"},
  {164.93032,"Ho","Holmium"},
  {167.259,"Er","Erbium"},
  {168.93421,"Tm","Thulium"},
  {173.04,"Yb","Ytterbium"},
  {174.967,"Lu","Lutetium"},
  {178.49,"Hf","Hafnium"},
  {180.9479,"Ta","Tantalum"},
  {183.84,"W","Tungsten"},
  {186.207,"Re","Rhenium"},
  {190.23,"Os","Osmium"},
  {192.217,"Ir","Iridium"},
  {195.078,"Pt","Platinum"},
  {196.96655,"Au","Gold"},
  {200.59,"Hg","Mercury"},
  {204.3833,"Tl","Thallium"},
  {207.2,"Pb","Lead"},
  {208.98038,"Bi","Bismuth"},
  {231.03588,"Pa","Protactinium"},
  {232.0381,"Th","Thorium"},
  {238.02891,"U","Uranium"}
};

static const int nelements = sizeof(amu)/sizeof(amu[0]);

/* look up element by name, case insensitive.  Return atomic number if
 * found, 0 if not found */
static int find_element(const char* abbrv) {
    for (int i=0; i<nelements; i++) {
        if (!strcasecmp(abbrv, amu[i].abbreviation)) {
            return i+1;
        }
    }
    return 0;
}

void adef_t::parse(const char *aref) {
    if ( isdigit(*aref) ) { res = *aref - '1'; ++aref; }
    else { res = 0; }
    if ( *aref == '-' ) { rel = -1; ++aref; }
    else if ( *aref == '+' ) { rel = 1; ++aref; }
    else if ( *aref == '#' ) { rel = 2; ++aref; }
    else { rel = 0; }
    name = aref;
}

static void print_msg(void *, const char *s) {
    fprintf(stdout, "%s\n", s);
}

static void debug_msg(const char *s) {
    if (getenv("MSYS_BUILDER_DEBUG")) {
        fprintf(stderr, "%s\n", s);
    }
}

void defs_t::import_charmm_topology(std::string const& path) {

    char *tok[TOKLEN];
    char sbuf[BUFLEN];
    char msgbuf[BUFLEN];
    int ntok;
    int itok;
    int first;
    int skip;
    int skipall;
    int stream;

    unsigned int utmp;

    FILE* file = fopen(path.c_str(), "r");
    if (!file) throw std::runtime_error("error opening topology file");
    boost::shared_ptr<FILE> ptr(file, fclose);
    
    static const bool all_caps = false;
    static void* v = NULL;

    first = 1;
    skip = 0;
    skipall = 0;
    stream = 0;

    resdef_t * res = NULL;

    while ( (ntok = charmm_get_tokens(tok,TOKLEN,sbuf,BUFLEN,file,all_caps)) ) {
        if ( ! tok[0][0] ) {
            print_msg(v,tok[1]);
            continue;
        }
        if ( skipall ) {
            print_msg (v, "skipping statements at end of file due to end or return statement");
            break;
        }
        if ( first ) {
            first = 0;
            if ( ! strncasecmp("READ",tok[0],4) && ! strncasecmp("RTF",tok[1],4) ) {
                print_msg (v, "reading topology from stream file");
                first = 1;
                stream = 1;
                continue;
            } else if ( ! strncasecmp("READ",tok[0],4) && ! strncasecmp("PARA",tok[1],4) ) {
                print_msg (v, "skipping parameters in stream file");
                skip = 1;
                continue;
            } else if ( ! strncasecmp("READ",tok[0],4) ) {
                print_msg (v, "skipping unknown section in stream file");
                skip = 1;
                continue;
            } else if ( ntok == 2 && sscanf(tok[0],"%u",&utmp) == 1
                    && sscanf(tok[1],"%u",&utmp) == 1 ) {
                sprintf(msgbuf,"Created by CHARMM version %s %s",tok[0],tok[1]);
                print_msg(v,msgbuf);
                continue;
            }
        }

        if ( skip ) {
            if ( ! strncasecmp("END",tok[0],4) ) {
                debug_msg("Recognized file end statement in skipped section.");
                skip = 0;
                first = 1;
            }
        }
        else if ( ! strncasecmp("END",tok[0],4) ) {
            debug_msg("Recognized file end statement.");
            if ( stream ) {
                stream = 0;
                first = 1;
            } else {
                skipall = 1;
            }
        }
        else if ( ! strncasecmp("RETURN",tok[0],4) ) {
            debug_msg("Recognized return statement.");
            skipall = 1;
        }
        else if ( ! strncasecmp("ACCE",tok[0],4) ) {
            debug_msg("Recognized acceptor statement.");
        }
        else if ( ! strncasecmp("DONO",tok[0],4) ) {
            debug_msg("Recognized donor statement.");
        }
        else if ( ! strncasecmp("BOND",tok[0],4) ||
                ! strncasecmp("DOUB",tok[0],4) ||
                ! strncasecmp("TRIP",tok[0],4) ) {
            debug_msg("Recognized bond statement.");
            if ( ntok < 3 || (ntok-1)%2 || !res) {
                print_msg(v,"ERROR!  Failed to parse bond statement.");
            } else {
                for ( itok = 1; itok < ntok; itok += 2 ) {
                    bond_t& bond = res->add_bond();
                    bond.def1.parse(tok[itok]);
                    bond.def2.parse(tok[itok+1]);
                }
            }
        }
        else if ( ! strncasecmp("ANGL",tok[0],4) ||
                ! strncasecmp("THET",tok[0],4) ) {
            debug_msg("Recognized angle statement.");
            if ( ntok < 4 || (ntok-1)%3 ) {
                print_msg(v,"ERROR!  Failed to parse angle statement.");
            } else {
#if 0
                for ( itok = 1; itok < ntok; itok += 3 ) {
                    s1 = parse_atom(tok[itok],&i1,&j1);
                    s2 = parse_atom(tok[itok+1],&i2,&j2);
                    s3 = parse_atom(tok[itok+2],&i3,&j3);
                    if ( topo_defs_angle(defs,0,0,s1,i1,j1,s2,i2,j2,s3,i3,j3) )
                        print_msg(v,"ERROR!  Failed to parse angle statement.");
                }
#endif
            }
        }
        else if ( ! strncasecmp("DIHE",tok[0],4) ) {
            debug_msg("Recognized dihedral statement.");
            if ( ntok < 5 || (ntok-1)%4 ) {
                print_msg(v,"ERROR!  Failed to parse dihedral statement.");
            } else {
#if 0
                for ( itok = 1; itok < ntok; itok += 4 ) {
                    s1 = parse_atom(tok[itok],&i1,&j1);
                    s2 = parse_atom(tok[itok+1],&i2,&j2);
                    s3 = parse_atom(tok[itok+2],&i3,&j3);
                    s4 = parse_atom(tok[itok+3],&i4,&j4);
                    if (topo_defs_dihedral(defs,0,0,s1,i1,j1,s2,i2,j2,s3,i3,j3,s4,i4,j4))
                        print_msg(v,"ERROR!  Failed to parse dihedral statement.");
                }
#endif
            }
        }
        else if ( ! strncasecmp("IMPH",tok[0],4) ||
                ! strncasecmp("IMPR",tok[0],4) ) {
            debug_msg("Recognized improper statement.");
            if ( ntok < 5 || (ntok-1)%4 ) {
                print_msg(v,"ERROR!  Failed to parse improper statement.");
            } else {
#if 0
                for ( itok = 1; itok < ntok; itok += 4 ) {
                    s1 = parse_atom(tok[itok],&i1,&j1);
                    s2 = parse_atom(tok[itok+1],&i2,&j2);
                    s3 = parse_atom(tok[itok+2],&i3,&j3);
                    s4 = parse_atom(tok[itok+3],&i4,&j4);
                    if (topo_defs_improper(defs,0,0,s1,i1,j1,s2,i2,j2,s3,i3,j3,s4,i4,j4))
                        print_msg(v,"ERROR!  Failed to parse improper statement.");
                }
#endif
            }
        }
        else if ( ! strncasecmp("CMAP",tok[0],4) ) {
            debug_msg("Recognized CMAP statement.");
            if ( ntok != 9 ) {
                print_msg(v,"ERROR!  Failed to parse CMAP statement.");
            } else {
#if 0
                const char* s[8]; int i[8], j[8];
                for ( itok = 0; itok < 8; ++itok ) {
                    s[itok] = parse_atom(tok[itok+1],&i[itok],&j[itok]);
                }
                if (topo_defs_cmap(defs,0,0,s,i,j))
                    print_msg(v,"ERROR!  Failed to parse CMAP statement.");
#endif
            }
        }
        else if ( ! strncasecmp("DECL",tok[0],4) ) {
            debug_msg("Recognized atom declaration statement.");
        }
        else if ( ! strncasecmp("ATOM",tok[0],4) ) {
            debug_msg("Recognized atom statement.");
            if ( ntok < 4 || !res ) {
                print_msg(v,"ERROR!  Failed to parse atom statement.");
            } else if ( ntok > 4 ) {
                print_msg(v,"ERROR!  Explicit exclusions or fluctuating charges not supported, atom ignored.");
            } else {
                atom_t& atom = res->add_atom(tok[1]);
                atom.type = tok[2];
                atom.charge = atof(tok[3]);
            }
        }
        else if ( ! strncasecmp("MASS",tok[0],4) ) {
            debug_msg("Recognized mass statement.");
            if ( ntok < 4 ) {
                print_msg(v,"ERROR!  Failed to parse mass statement.");
            } else {
                type_t& type = add_type(tok[2]);
                const char* abbrv = ntok>4?tok[4]:"";
                type.anum = find_element(abbrv);
                type.mass = atof(tok[3]);
                type.id = atoi(tok[1]);
            }
        }
        else if ( ! strncasecmp("AUTO",tok[0],4) ) {
            debug_msg("Recognized autogenerate statement.");
#if 0
            for ( itok = 1; itok < ntok; itok += 1 ) {
                if ( ! strncasecmp("ANGL",tok[itok],4) ) {
                    topo_defs_auto_angles(defs,1);
                } else if ( ! strncasecmp("DIHE",tok[itok],4) ) {
                    topo_defs_auto_dihedrals(defs,1);
                } else {
                    print_msg(v,"ERROR!  Failed to parse autogenerate statement.");
                }
            }
#endif
        }
        else if ( ! strncasecmp("DEFA",tok[0],4) ) {
            debug_msg("Recognized default patch statement.");
            if ( ntok < 3 || (ntok-1)%2 ) {
                print_msg(v,"ERROR!  Failed to parse default patching statement.");
            } else {
                for ( itok = 1; itok < ntok; itok += 2 ) {
                    if ( ! strncasecmp("FIRS",tok[itok],4) ) {
                        this->pfirst = tok[itok+1];
                    } else if ( ! strncasecmp("LAST",tok[itok],4) ) {
                        this->plast = tok[itok+1];
                    } else {
                        print_msg(v,"ERROR!  Failed to parse default patching statement.");
                    }
                }
            }
        }
        else if ( ! strncasecmp("BILD",tok[0],4) ||
                ! strncasecmp("IC",tok[0],4) ) {
            debug_msg("Recognized internal coordinate statement.");
            if ( ntok < 10 || !res) {
                print_msg(v,"ERROR!  Failed to parse internal coordinate statement.");
            } else {
                conf_t& conf = res->add_conf();
                conf.improper = *tok[3]=='*';
                conf.def1.parse(tok[1]);
                conf.def2.parse(tok[2]);
                conf.def3.parse(tok[3] + conf.improper);
                conf.def4.parse(tok[4]);
                conf.r1 = atof(tok[5]);
                conf.a1 = atof(tok[6]) * M_PI / 180;
                conf.phi = atof(tok[7]) * M_PI / 180;
                conf.a2 = atof(tok[8]) * M_PI / 180;
                conf.r2 = atof(tok[9]);
            }
        }
        else if ( ! strncasecmp("DELE",tok[0],4) ) {
            debug_msg("Recognized delete statement.");
            if ( ntok < 2 || !res ) {
                print_msg(v,"ERROR!  Failed to parse delete statement.");
            } else {
                if ( ! strncasecmp("ATOM",tok[1],4) ) {
                    if ( ntok < 3 ) {
                        print_msg(v,"ERROR!  Failed to parse delete atom statement.");
                    } else {
                        adef_t& atom = res->add_delatom();
                        atom.parse(tok[2]);
                    }
                } else if ( ! strncasecmp("ACCE",tok[1],4) ) {
                    ;
                } else if ( ! strncasecmp("DONO",tok[1],4) ) {
                    ;
                } else if ( ! strncasecmp("BOND",tok[1],4) ||
                        ! strncasecmp("DOUB",tok[1],4) ||
                        ! strncasecmp("TRIP",tok[1],4) ) {
                    if ( ntok < 4 || (ntok-2)%2 ) {
                        print_msg(v,"ERROR!  Failed to parse delete bond statement.");
                    } else {
                        for ( itok = 2; itok < ntok; itok += 2 ) {
                            bond_t& bond = res->add_delbond();
                            bond.def1.parse(tok[itok]);
                            bond.def2.parse(tok[itok+1]);
                        }
                    }
                } else if ( ! strncasecmp("ANGL",tok[1],4) ||
                        ! strncasecmp("THET",tok[1],4) ) {
                    if ( ntok < 5 || (ntok-2)%3 ) {
                        print_msg(v,"ERROR!  Failed to parse delete angle statement.");
                    } else {
#if 0
                        for ( itok = 2; itok < ntok; itok += 3 ) {
                            s1 = parse_atom(tok[itok],&i1,&j1);
                            s2 = parse_atom(tok[itok+1],&i2,&j2);
                            s3 = parse_atom(tok[itok+2],&i3,&j3);
                            if ( topo_defs_angle(defs,0,1,s1,i1,j1,s2,i2,j2,s3,i3,j3) )
                                print_msg(v,"ERROR!  Failed to parse delete angle statement.");
                        }
#endif
                    }
                } else if ( ! strncasecmp("DIHE",tok[1],4) ) {
                    if ( ntok < 6 || (ntok-2)%4 ) {
                        print_msg(v,"ERROR!  Failed to parse delete dihedral statement.");
                    } else {
#if 0
                        for ( itok = 2; itok < ntok; itok += 4 ) {
                            s1 = parse_atom(tok[itok],&i1,&j1);
                            s2 = parse_atom(tok[itok+1],&i2,&j2);
                            s3 = parse_atom(tok[itok+2],&i3,&j3);
                            s4 = parse_atom(tok[itok+3],&i4,&j4);
                            if (topo_defs_dihedral(defs,0,1,s1,i1,j1,s2,i2,j2,s3,i3,j3,s4,i4,j4))
                                print_msg(v,"ERROR!  Failed to parse delete dihedral statement.");
                        }
#endif
                    }
                } else if ( ! strncasecmp("IMPH",tok[1],4) ||
                        ! strncasecmp("IMPR",tok[1],4) ) {
                    if ( ntok < 6 || (ntok-2)%4 ) {
                        print_msg(v,"ERROR!  Failed to parse delete improper statement.");
                    } else {
#if 0
                        for ( itok = 2; itok < ntok; itok += 4 ) {
                            s1 = parse_atom(tok[itok],&i1,&j1);
                            s2 = parse_atom(tok[itok+1],&i2,&j2);
                            s3 = parse_atom(tok[itok+2],&i3,&j3);
                            s4 = parse_atom(tok[itok+3],&i4,&j4);
                            if (topo_defs_improper(defs,0,1,s1,i1,j1,s2,i2,j2,s3,i3,j3,s4,i4,j4))
                                print_msg(v,"ERROR!  Failed to parse delete improper statement.");
                        }
#endif
                    }
                } else if ( ! strncasecmp("BILD",tok[1],4) ||
                        ! strncasecmp("IC",tok[1],4) ) {
                    if ( ntok < 6 ) {
                        print_msg(v,"ERROR!  Failed to parse delete internal coordinate statement.");
                    } else {
                        conf_t& conf = res->add_conf();
                        conf.improper = *tok[4]=='*';
                        conf.def1.parse(tok[2]);
                        conf.def2.parse(tok[3]);
                        conf.def3.parse(tok[4] + conf.improper);
                        conf.def4.parse(tok[5]);
                    }
                } else {
                    print_msg(v,"ERROR!  Failed to parse delete statement.");
                }
            }
        }
        else if ( ! strncasecmp("GROU",tok[0],4) ) {
            debug_msg("Recognized group statement.");
        }
        else if ( ! strncasecmp("PATC",tok[0],4) ) {
            debug_msg("Recognized patching statement.");
            if ( ntok < 3 || (ntok-1)%2 ) {
                print_msg(v,"ERROR!  Failed to parse patching statement.");
            } else {
                for ( itok = 1; itok < ntok; itok += 2 ) {
                    if ( ! strncasecmp("FIRS",tok[itok],4) ) {
                        res->pfirst = tok[itok+1];
                    } else if ( ! strncasecmp("LAST",tok[itok],4) ) {
                        res->plast = tok[itok+1];
                    } else {
                        print_msg(v,"ERROR!  Failed to parse patching statement.");
                    }
                }
            }
        }
        else if ( ! strncasecmp("RESI",tok[0],4) ) {
            debug_msg("Recognized residue statement.");
            if ( ntok < 2 ) {
                print_msg(v,"ERROR!  Failed to parse residue statement.");
            } else {
                res = &add_resdef(tok[1]);
            }
        }
        else if ( ! strncasecmp("PRES",tok[0],4) ) {
            debug_msg("Recognized patch residue statement.");
            if ( ntok < 2 ) {
                print_msg(v,"ERROR!  Failed to parse patch residue statement.");
            } else {
                res = &add_resdef(tok[1]);
                res->patch = true;
            }
        }
        else {
            sprintf(msgbuf,"ERROR!  FAILED TO RECOGNIZE %s",tok[0]);
            print_msg(v,msgbuf);
        }

    }

}

void resdef_t::patch_topology(resdef_t& topo) const {
    for (unsigned i=0; i<delatoms.size(); i++) {
        if (0==topo.atoms.erase(delatoms[i].name)) {
            printf("ERROR, failed to erase atom %s\n", delatoms[i].name.c_str());
        }
        /* delete bonds to-from that atom */
        unsigned j=topo.bonds.size();
        while (j) {
            --j;
            if (topo.bonds[j].def1.name==delatoms[i].name ||
                topo.bonds[j].def2.name==delatoms[i].name) {
                topo.bonds.erase(topo.bonds.begin()+j);
            }
        }
        /* delete confs to-from that atom */
        j=topo.confs.size();
        while (j) {
            --j;
            if (topo.confs[j].def1.name==delatoms[i].name ||
                topo.confs[j].def2.name==delatoms[i].name ||
                topo.confs[j].def3.name==delatoms[i].name ||
                topo.confs[j].def4.name==delatoms[i].name) {
                topo.confs.erase(topo.confs.begin()+j);
            }
        }
    }

    /* add structure from patch */
    topo.atoms.insert(atoms.begin(), atoms.end());
    topo.bonds.insert(topo.bonds.begin(), bonds.begin(), bonds.end());
    topo.confs.insert(topo.confs.begin(), confs.begin(), confs.end());
}
