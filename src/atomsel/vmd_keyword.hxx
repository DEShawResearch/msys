/* @COPYRIGHT@ */

#ifndef desres_msys_atomsel_vmd_keyword_hxx
#define desres_msys_atomsel_vmd_keyword_hxx

#include "keyword.hxx"
#include "../system.hxx"

namespace desres { namespace msys { namespace atomsel {

KeywordPtr keyword_fragment( SystemPtr sys );
KeywordPtr keyword_water( SystemPtr sys );
KeywordPtr keyword_backbone( SystemPtr sys );
KeywordPtr keyword_protein( SystemPtr sys );
KeywordPtr keyword_nucleic( SystemPtr sys );

}}}

#endif
