/* @COPYRIGHT@ */

#ifndef desres_msys_atomsel_msys_keyword_hxx
#define desres_msys_atomsel_msys_keyword_hxx

#include "keyword.hxx"
#include "../system.hxx"

namespace desres { namespace msys { namespace atomsel {


/* A selection of all the atoms in the system */
Selection full_selection(SystemPtr sys);

/* all true */
KeywordPtr keyword_all( SystemPtr ent );

/* all false */
KeywordPtr keyword_none( SystemPtr ent );


KeywordPtr keyword_alchemical( SystemPtr ent );
KeywordPtr keyword_anum( SystemPtr ent );
KeywordPtr keyword_chain( SystemPtr ent );
KeywordPtr keyword_element( SystemPtr ent );
KeywordPtr keyword_hydrogen( SystemPtr ent );
KeywordPtr keyword_index( SystemPtr ent );
KeywordPtr keyword_mass( SystemPtr ent );
KeywordPtr keyword_charge( SystemPtr ent );
KeywordPtr keyword_name( SystemPtr ent );
KeywordPtr keyword_resid( SystemPtr ent );
KeywordPtr keyword_residue( SystemPtr ent );
KeywordPtr keyword_resname( SystemPtr ent );
KeywordPtr keyword_numbonds( SystemPtr ent );
KeywordPtr keyword_residue( SystemPtr ent );
KeywordPtr keyword_fragid( SystemPtr ent );
KeywordPtr keyword_x( SystemPtr ent );
KeywordPtr keyword_y( SystemPtr ent );
KeywordPtr keyword_z( SystemPtr ent );
KeywordPtr keyword_vx( SystemPtr ent );
KeywordPtr keyword_vy( SystemPtr ent );
KeywordPtr keyword_vz( SystemPtr ent );
KeywordPtr keyword_water( SystemPtr sys );
KeywordPtr keyword_backbone( SystemPtr sys );
KeywordPtr keyword_protein( SystemPtr sys );
KeywordPtr keyword_nucleic( SystemPtr sys );

KeywordPtr keyword_atomprop( SystemPtr ent, String const& prop );

}}}

#endif
