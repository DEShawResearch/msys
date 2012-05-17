/***************************************************************************
 *cr
 *cr            (C) Copyright 1995-2009 The Board of Trustees of the
 *cr                        University of Illinois
 *cr                         All Rights Reserved
 *cr
 ***************************************************************************/

/***************************************************************************
 * RCS INFORMATION:
 *
 *      $RCSfile: readpdb.h,v $
 *      $Author$       $Locker:  $             $State: Exp $
 *      $Revision$       $Date$
 *
 ***************************************************************************/

#ifndef READ_PDB_H
#define READ_PDB_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define PDB_RECORD_LENGTH   80   /* actual record size */
#define PDB_BUFFER_LENGTH   83   /* size need to buffer + CR, LF, and NUL */

#define VMDUSECONECTRECORDS 1

#ifdef __cplusplus
extern "C" {
#endif

/*  record type defines */
enum {
  PDB_HEADER, PDB_REMARK, PDB_ATOM, PDB_CONECT, PDB_UNKNOWN, PDB_END, PDB_EOF, PDB_CRYST1
};

/* read the next record from the specified pdb file, and put the string found
   in the given string pointer (the caller must provide adequate (81 chars)
   buffer space); return the type of record found
*/
int desres_msys_read_pdb_record(FILE *f, char *retStr);

/* Extract the alpha/beta/gamma a/b/c unit cell info from a CRYST1 record */
void desres_msys_get_pdb_cryst1(const char *record, 
                           float *alpha, float *beta, float *gamma, 
                           float *a, float *b, float *c);

/* Extract the x,y,z coords, occupancy, and beta from an ATOM record */
void desres_msys_get_pdb_coordinates(const char *record, 
                         float *x, float *y, float *z,
                         float *occup, float *beta);

void desres_msys_get_pdb_header(const char *record, char *pdbcode, char *date,
                           char *classification);

void desres_msys_get_pdb_conect(const char *record, int natoms, int *idxmap,
                           int *maxbnum, int *nbonds, int **from, int **to);

/* ATOM field format according to PDB standard v2.2
  COLUMNS        DATA TYPE       FIELD         DEFINITION
---------------------------------------------------------------------------------
 1 -  6        Record name     "ATOM  "
 7 - 11        Integer         serial        Atom serial number.
13 - 16        Atom            name          Atom name.
17             Character       altLoc        Alternate location indicator.
18 - 20        Residue name    resName       Residue name.
22             Character       chainID       Chain identifier.
23 - 26        Integer         resSeq        Residue sequence number.
27             AChar           iCode         Code for insertion of residues.
31 - 38        Real(8.3)       x             Orthogonal coordinates for X in Angstroms.
39 - 46        Real(8.3)       y             Orthogonal coordinates for Y in Angstroms.
47 - 54        Real(8.3)       z             Orthogonal coordinates for Z in Angstroms.
55 - 60        Real(6.2)       occupancy     Occupancy.
61 - 66        Real(6.2)       tempFactor    Temperature factor.
73 - 76        LString(4)      segID         Segment identifier, left-justified.
77 - 78        LString(2)      element       Element symbol, right-justified.
79 - 80        LString(2)      charge        Charge on the atom.
 */

/* Break a pdb ATOM record into its fields.  The user must provide the
   necessary space to store the atom name, residue name, and segment name.
   Character strings will be null-terminated.
*/
void desres_msys_get_pdb_fields(const char *record, int reclength, int *serial,
                    char *name, char *resname, char *chain, 
                    char *segname, char *resid, char *insertion, 
                    char *altloc, char *elementsymbol,
                    float *x, float *y, float *z, 
                    float *occup, float *beta);


#ifdef __cplusplus
}
#endif

#endif
