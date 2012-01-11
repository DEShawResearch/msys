#ifndef desres_msys_builder_builder_hxx
#define desres_msys_builder_builder_hxx

#include "defs.hxx"
#include "../../src/system.hxx"

namespace desres { namespace msys { namespace builder {

    void build( defs_t const& defs, SystemPtr mol, Id chain,
                std::string pfirst="", std::string plast="");

    void patch( defs_t const& defs, std::string const& pres, 
                SystemPtr mol, IdList const& residues );


}}}

#endif
