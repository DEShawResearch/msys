#include "../psf.hxx"
#include "../schema.hxx"
#include "../elements.hxx"
#include "../term_table.hxx"
#include <stdio.h>

//using namespace desres::msys;


#define PSF_RECORD_LENGTH 256  /* extended to handle Charmm CMAP/CHEQ/DRUDE */

/* Formatted reads:
 *
 * copy at most 'len' characters from source to target. 
 *
 * leading white space is skipped over but counts towards 'len'.
 * the copy stops at first whitspace or a '\0'.
 * unlike strncpy(3) the result will always be \0 terminated.
 *
 * intended for copying (short) strings from formatted fortran  
 * i/o files that must not contain whitespace (e.g. residue names, 
 * atom name/types etc. in .pdb, .psf and alike.). 
 */
static void strnwscpy(char *target, const char *source, const int len) {
  int i, c;

  for (i=0, c=0; i<len; ++i) {
    if (*source == '\0' || (c > 0 && *source == ' ')) {
      break;
    }

    if (*source == ' ') { 
      source++;
    } else {
      *target++ = *source++;
      c++;
    }
  }
  *target = '\0';
}

/* Formatted reads:
 *
 * copy at most 'maxlen' characters from source to target allowing overflow.
 *
 * leading white space up to 'len' is skipped over but counts towards 'maxlen'.
 * the copy stops at first whitspace or a '\0'.
 * unlike strncpy(3) the result will always be \0 terminated.
 *
 * intended for copying (short) strings from formatted fortran
 * i/o files that must not contain whitespace (e.g. residue names,
 * atom name/types etc. in .pdb, .psf and alike.).
 *
 * returns number of bytes of overflow.
 */
static int strnwscpy_shift(char *target, const char *source,
                           const int len, const int maxlen) {
  int i, c;

  for (i=0, c=0; i<maxlen; ++i) {
    if (*source == '\0' || (c > 0 && *source == ' ') || (c == 0 && i == len)) {
      break;
    }

    if (*source == ' ') {
      source++;
    } else {
      *target++ = *source++;
      c++;
    }
  }
  *target = '\0';
  return ( i > len ? i - len : 0 );
}

/* atoi() replacement
 *
 * reads int with field width fw handling various overflow cases to
 * support both " %7d %7d" and "%8d%8d" writers up to 100M atoms.
 *
 */

static int atoifw(char **ptr, int fw) {
  char *op = *ptr;
  int ival = 0;
  int iws = 0;
  char tmpc;

  sscanf(op, "%d%n", &ival, &iws);
  if ( iws == fw ) { /* "12345678 123..." or " 1234567 123..." */
    *ptr += iws;
  } else if ( iws < fw ) { /* left justified? */
    while ( iws < fw && op[iws] == ' ' ) ++iws;
    *ptr += iws;
  } else if ( iws < 2*fw ) { /* " 12345678 123..." */
    *ptr += iws;
  } else { /* " 123456712345678" or "1234567812345678" */
    tmpc = op[fw];  op[fw] = '\0';
    ival = atoi(op);
    op[fw] = tmpc;
    *ptr += fw;
  }
  return ival;
}


/* Read in the next atom info into the given storage areas; this assumes
   that file has already been moved to the beginning of the atom records.
   Returns the serial number of the atom. If there is an error, returns -1.*/
static int get_psf_atom(FILE *f, char *name, char *atype, char *resname,
                        char *segname, int *resid, char *insertion, double *q, double *m, 
                        int namdfmt, int charmmext, int charmmdrude) {
  char inbuf[PSF_RECORD_LENGTH+2];
  int num;

  if (inbuf != fgets(inbuf, PSF_RECORD_LENGTH+1, f)) {
    return(-1); /* failed to read in an atom */
  }

  if (strlen(inbuf) < 50) {
    fprintf(stderr, "Line too short in psf file: \n%s\n", inbuf);
    return -1;
  }

  num = atoi(inbuf); /* atom index */

  if (namdfmt == 1) {
    int cnt, rcnt;
    char residstr[8], trash;
    cnt = sscanf(inbuf, "%d %7s %7s %7s %7s %7s %lf %lf",
                 &num, segname, residstr, resname, name, atype, q, m);
    insertion[0] = ' ';  insertion[1] = '\0';
    rcnt = sscanf(residstr, "%d%c%c", resid, insertion, &trash);
    if (cnt != 8 || rcnt < 1 || rcnt > 2) {
      printf("psfplugin) Failed to parse atom line in NAMD PSF file:\n");
      printf("psfplugin)   '%s'\n", inbuf);
      return -1;
    }
  } else if (charmmdrude == 1 || charmmext == 1) {
    int xplorshift;
    /* CHARMM PSF format is (if DRUDE or (?) CHEQ are enabled):
     *  '(I10,1X,A8,1X,A8,1X,A8,1X,A8,1X,I4,1X,2G14.6,I8,2G14.6)'
     */
    if ( inbuf[10] != ' ' ||
         inbuf[19] != ' ' ||
         inbuf[28] != ' ' ||
         inbuf[37] != ' ' ||
         inbuf[46] != ' ' ) {
      printf("psfplugin) Failed to parse atom line in PSF file:\n");
      printf("psfplugin)   '%s'\n", inbuf);
      return -1;
    }

    strnwscpy(segname, inbuf+11, 7);
    strnwscpy(resname, inbuf+29, 7);
    strnwscpy(name, inbuf+38, 7);

    xplorshift = 0;
    strnwscpy(atype, inbuf+47, 4);
    if ( ! isdigit(atype[0]) ) {
      strnwscpy(atype, inbuf+47, 6);
      xplorshift = 2;
    }

    if ( inbuf[51+xplorshift] != ' ' ) {
      printf("psfplugin) Failed to parse atom line in PSF file:\n");
      printf("psfplugin)   '%s'\n", inbuf);
      return -1;
    }
    
    insertion[0] = ' ';  insertion[1] = '\0';
    sscanf(inbuf+20, "%d%c", resid, insertion);
    *q = atof(inbuf+52+xplorshift);
    *m = atof(inbuf+66+xplorshift);
    // data we don't currently read:
    // if (charmmdrude == 1) {
    //   *imove = atoi(inbuf+80+xplorshift);
    //   *alphadp = atof(inbuf+88+xplorshift);
    //   *tholei = atof(inbuf+102+xplorshift);
    // }
  } else {
    /* CHARMM PSF format is 
     *  '(I8,1X,A4,1X,A4,1X,A4,1X,A4,1X,I4,1X,2G14.6,I8)'
     */
    const char *rdbuf = inbuf;
    char intbuf[16];

    intbuf[0] = '\0';
    rdbuf += strnwscpy_shift(intbuf, rdbuf, 8, 10);
    if ( rdbuf[8] != ' ' ) {
      printf("psfplugin) Failed to parse atom index in PSF file:\n");
      printf("psfplugin)   '%s'\n", inbuf);
      return -1;
    }
    rdbuf += strnwscpy_shift(segname, rdbuf+9, 4, 7);
    if ( rdbuf[13] != ' ' ) {
      printf("psfplugin) Failed to parse segname in PSF file:\n");
      printf("psfplugin)   '%s'\n", inbuf);
      return -1;
    }
    intbuf[0] = '\0';
    rdbuf += strnwscpy_shift(intbuf, rdbuf+14, 4, 8);
    insertion[0] = ' ';  insertion[1] = '\0';
    sscanf(intbuf, "%d%c", resid, insertion);
    if ( rdbuf[18] != ' ' ) {
      printf("psfplugin) Failed to parse resid in PSF file:\n");
      printf("psfplugin)   '%s'\n", inbuf);
      return -1;
    }
    rdbuf += strnwscpy_shift(resname, rdbuf+19, 4, 7);
    if ( rdbuf[23] != ' ' ) {
      printf("psfplugin) Failed to parse resname in PSF file:\n");
      printf("psfplugin)   '%s'\n", inbuf);
      return -1;
    }
    rdbuf += strnwscpy_shift(name, rdbuf+24, 4, 7);
    if ( rdbuf[28] != ' ' ) {
      printf("psfplugin) Failed to parse atom name in PSF file:\n");
      printf("psfplugin)   '%s'\n", inbuf);
      return -1;
    }
    rdbuf += strnwscpy_shift(atype, rdbuf+29, 4, 7);
    if ( rdbuf[33] != ' ' ) {
      printf("psfplugin) Failed to parse atom type in PSF file:\n");
      printf("psfplugin)   '%s'\n", inbuf);
      return -1;
    }
    *q = (float) atof(rdbuf+34);
    *m = (float) atof(rdbuf+48);
  }

#if 0
  /* if this is a Charmm31 PSF file, there may be two extra */
  /* columns containing polarizable force field data.       */
  if (psf->charmmcheq) {
    /* do something to read in these columns here */
  }
#endif

  return num;
}


/*
 * Read in the beginning of the bond/angle/dihed/etc information,
 * but don't read in the data itself.  Returns the number of the record type
 * for the molecule.  If error, returns (-1). 
 */
static int psf_start_block(FILE *file, const char *blockname) {
  char inbuf[PSF_RECORD_LENGTH+2];
  int nrec = -1;
  
  /* check if we had a parse error earlier, which is indicated
     by the file descriptor set to NULL */
  if (!file)
    return -1;

  /* keep reading the next line until a line with blockname appears */
  do {
    if(inbuf != fgets(inbuf, PSF_RECORD_LENGTH+1, file)) {
      /* EOF encountered with no blockname line found ==> error, return (-1) */
      return -1;
    }
    if(strlen(inbuf) > 0 && strstr(inbuf, blockname))
      nrec = atoi(inbuf);
  } while (nrec == -1);

  return nrec;
}


/* Read in the bond info into the given integer arrays, one for 'from' and
   one for 'to' atoms; remember that .psf files use 1-based indices,
   not 0-based.  Returns 1 if all nbond bonds found; 0 otherwise.  */
static int psf_get_bonds(FILE *f, int nbond, int fromAtom[], int toAtom[], int charmmext, int namdfmt) {
  char *bondptr=NULL;
  int fw = charmmext ? 10 : 8;
  char inbuf[PSF_RECORD_LENGTH+2];
  int i=0;
  size_t minlinesize;
  int rc=0;

  while (i < nbond) {
    if (namdfmt) {
      // NAMD assumes a space-delimited variant of the PSF file format
      int cnt = fscanf(f, "%d %d", &fromAtom[i], &toAtom[i]);
      if (cnt < 2) {
        fprintf(stderr, "Bonds line too short in NAMD psf file.\n");
        break;
      }
    } else {
      if ((i % 4) == 0) {
        /* must read next line */
        if (!fgets(inbuf, PSF_RECORD_LENGTH+2, f)) {
          /* early EOF encountered */
          break;
        }

        /* Check that there is enough space in the line we are about to read */
        if (nbond-i >= 4) {
          minlinesize = 2*fw*4; 
        } else {
          minlinesize = 2*fw*(nbond-i); 
        }

        if (strlen(inbuf) < minlinesize) {
          fprintf(stderr, "Bonds line too short in psf file: \n%s\n", inbuf);
          break;
        }
        bondptr = inbuf;
      }

      if ((fromAtom[i] = atoifw(&bondptr,fw)) < 1) {
        printf("psfplugin) ERROR: Bond %d references atom with index < 1!\n", i);
        rc=-1;
        break;
      }
  
      if ((toAtom[i] = atoifw(&bondptr,fw)) < 1) {
        printf("psfplugin) ERROR: Bond %d references atom with index < 1!\n", i);
        rc=-1;
        break;
      }
    }

    i++;
  }

  if (rc == -1) {
    printf("psfplugin) ERROR: skipping bond info due to bad atom indices\n");
  } else if (i != nbond) {
    printf("psfplugin) ERROR: unable to read the specified number of bonds!\n");
    printf("psfplugin) Expected %d bonds but only read %d\n", nbond, i);
  }

  return (i == nbond);
}

static int psf_get_angles(FILE *f, int n, int *angles, int charmmext) {
    char inbuf[PSF_RECORD_LENGTH+2];
    char *bondptr = NULL;
    int fw = charmmext ? 10 : 8;
    int i=0;
    while (i<n) {
        if((i % 3) == 0) {
            /* must read next line */
            if(!fgets(inbuf,PSF_RECORD_LENGTH+2,f)) {
                /* early EOF encountered */
                break;
            }
            bondptr = inbuf;
        }
        if((angles[3*i] = atoifw(&bondptr,fw)) < 1)
            break;
        if((angles[3*i+1] = atoifw(&bondptr,fw)) < 1)
            break;
        if((angles[3*i+2] = atoifw(&bondptr,fw)) < 1)
            break;
        i++;
    }

    return (i != n);
}

static int psf_get_dihedrals_impropers(FILE *f, int n, int *dihedrals, int charmmext) {
    char inbuf[PSF_RECORD_LENGTH+2];
    char *bondptr = NULL;
    int fw = charmmext ? 10 : 8;
    int i=0;
    while (i<n) {
        if((i % 2) == 0) {
            /* must read next line */
            if(!fgets(inbuf,PSF_RECORD_LENGTH+2,f)) {
                /* early EOF encountered */
                break;
            }
            bondptr = inbuf;
        }
        if((dihedrals[4*i] = atoifw(&bondptr,fw)) < 1)
            break;
        if((dihedrals[4*i+1] = atoifw(&bondptr,fw)) < 1)
            break;
        if((dihedrals[4*i+2] = atoifw(&bondptr,fw)) < 1)
            break;
        if((dihedrals[4*i+3] = atoifw(&bondptr,fw)) < 1)
            break;
        i++;
    }

    return (i != n);
}

desres::msys::SystemPtr desres::msys::ImportPSF( std::string const& path ) {

    FILE *fp;
    char inbuf[PSF_RECORD_LENGTH*8+2];
    const char *progname = "Charmm";

    /* Open the .psf file and skip past the remarks to the first data section.
     * Returns the file pointer, or NULL if error.  Also puts the number of
     * atoms in the molecule into the given integer.  
     */

    if ((fp = fopen(path.c_str(), "r")) == NULL) {
        MSYS_FAIL("Couldn't open psf file at " << path);
    }
    boost::shared_ptr<FILE> fp_closer(fp, fclose);

    bool namdfmt = 0;   /* off unless we discover otherwise */
    //bool charmmfmt = 0; /* off unless we discover otherwise */
    bool charmmext = 0; /* off unless we discover otherwise */
    bool charmmcheq = 0;
    bool charmmcmap = 0;
    bool charmmdrude = 0;
    int natoms = -1;

    /* read lines until a line with NATOM and without REMARKS appears    */
    do {
        /* be prepared for long lines from CNS remarks */
        if (inbuf != fgets(inbuf, PSF_RECORD_LENGTH*8+1, fp)) {
            /* EOF encountered with no NATOM line found ==> error */
            MSYS_FAIL("unexpected EOF in psf file");
        }

        if (strlen(inbuf) > 0) {
            if (!strstr(inbuf, "REMARKS")) {
                if (strstr(inbuf, "PSF")) {
                    if (strstr(inbuf, "NAMD")) {
                        namdfmt = 1;      
                    }
                    if (strstr(inbuf, "EXT")) {
                        //charmmfmt = 1; 
                        charmmext = 1;      
                    }
                    if (strstr(inbuf, "CHEQ")) {
                        //charmmfmt = 1; 
                        charmmcheq = 1;      
                    }
                    if (strstr(inbuf, "CMAP")) {
                        //charmmfmt = 1; 
                        charmmcmap = 1;      
                    }
                    if (strstr(inbuf, "DRUDE")) {
                        //charmmfmt = 1; 
                        charmmdrude = 1;      
                    }
                } else if (strstr(inbuf, "NATOM")) {
                    natoms = atoi(inbuf);
                }
            } 
        }
    } while (natoms == -1);

    if (namdfmt) {
        progname = "NAMD";
    } else {
        progname = "Charmm";
    }
    if (charmmcheq || charmmcmap) {
        printf("psfplugin) Detected a %s PSF file\n", progname);
    }
    if (charmmext) {
        printf("psfplugin) Detected a %s PSF EXTEnded file\n", progname);
    }
    if (charmmdrude) {
        printf("psfplugin) Detected a %s Drude polarizable force field file\n", progname);
        printf("psfplugin) WARNING: Support for Drude FF is currently experimental\n");
    }

    SystemPtr mol = System::create();
    SystemImporter imp(mol);

    typedef std::map<std::string, Id> TypeMap;
    TypeMap types;
    TermTablePtr nb = AddNonbonded(mol, "vdw_12_6", "arithmetic/geometric");
    ParamTablePtr params = nb->params();
    Id col = params->addProp("type", StringType);

    /* read atoms */
    for (int i=0; i<natoms; i++) {
        char name[8];
        char type[8];
        char resname[8];
        char segid[8];
        char chain[2];
        int resid;
        char insertion[4];
        double charge;
        double mass;

        if (get_psf_atom(fp, name, type, resname, segid, 
                    &resid, insertion, &charge, &mass,
                    namdfmt, charmmext, charmmdrude) < 0) {
            MSYS_FAIL("Error reading psf atom " << i+1);
        }
        chain[0]=segid[0];
        chain[1]='\0';
        Id id = imp.addAtom(chain, segid, resid, resname, name, insertion);
        atom_t& atm = mol->atomFAST(id);
        atm.mass = mass;
        atm.charge = charge;
        if(atm.mass>1.0) atm.atomic_number=GuessAtomicNumber(atm.mass);
        if (types.find(type)==types.end()) {
            Id param = params->addParam();
            types[type] = param;
            params->value(param, col) = type;
        } 
        nb->addTerm(IdList(1,id), types[type]);
    }

    /* read bonds */
    int nbonds = psf_start_block(fp, "NBOND"); /* get bond count */

    if (nbonds > 0) {
        TermTablePtr stretch = AddTable(mol, "stretch_harm");
        Id param = stretch->params()->addParam();
        IdList ids(2);
        std::vector<int> from(nbonds), to(nbonds);
        if (!psf_get_bonds(fp, nbonds, &from[0], &to[0],
                    charmmext, namdfmt)) {
            MSYS_FAIL("Error reading bonds");
        }
        for (int i=0; i<nbonds; i++) {
            mol->addBond(from[i]-1, to[i]-1);
            ids[0]=from[i]-1;
            ids[1]=to[i]-1;
            stretch->addTerm(ids, param);
        }
    }

    /* read angles */
    int nangles = psf_start_block(fp, "NTHETA"); /* get angle count */
    if (nangles > 0) {
        TermTablePtr table = AddTable(mol, "angle_harm");
        Id param = table->params()->addParam();
        IdList tmp(3);
        std::vector<int> ids(3*nangles);
        if (psf_get_angles(fp, nangles, &ids[0], charmmext)) {
            MSYS_FAIL("Error reading angles");
        }
        for (int i=0; i<nangles; i++) {
            for (int j=0; j<3; j++) tmp[j] = ids[3*i+j]-1;
            table->addTerm(tmp, param);
        }
    }

    /* read dihedrals */
    int ndihed = psf_start_block(fp, "NPHI");
    if (ndihed > 0) {
        TermTablePtr table = AddTable(mol, "dihedral_trig");
        Id param = table->params()->addParam();
        IdList tmp(4);
        std::vector<int> ids(4*ndihed);
        if (psf_get_dihedrals_impropers(fp, ndihed, &ids[0], charmmext)) {
            MSYS_FAIL("Error reading dihedrals");
        }
        for (int i=0; i<ndihed; i++) {
            for (int j=0; j<4; j++) tmp[j] = ids[4*i+j]-1;
            table->addTerm(tmp, param);
        }
    }

    /* read impropers */
    int nimp = psf_start_block(fp, "NIMPHI");
    if (nimp > 0) {
        TermTablePtr table = AddTable(mol, "improper_harm");
        Id param = table->params()->addParam();
        IdList tmp(4);
        std::vector<int> ids(4*nimp);
        psf_get_dihedrals_impropers(fp, nimp, &ids[0], charmmext);
        for (int i=0; i<nimp; i++) {
            for (int j=0; j<4; j++) tmp[j] = ids[4*i+j]-1;
            table->addTerm(tmp, param);
        }
    }

    /* read cross-terms (cmap) */
    int ncmap = psf_start_block(fp, "NCRTERM");
    if (ncmap> 0) {
        TermTablePtr table = AddTable(mol, "torsiontorsion_cmap");
        Id param = table->params()->addParam();
        IdList tmp(8);
        std::vector<int> ids(8*ncmap);
        psf_get_dihedrals_impropers(fp, 2*ncmap, &ids[0], charmmext);
        for (int i=0; i<ncmap; i++) {
            for (int j=0; j<8; j++) tmp[j] = ids[8*i+j]-1;
            table->addTerm(tmp, param);
        }
    }

    return mol;
}

