#ifndef msys_atomsel_k_nearest_hxx
#define msys_atomsel_k_nearest_hxx

#include "predicate.hxx"
#include "../system.hxx"

namespace desres { namespace msys { namespace atomsel {

  /*! true for the k nearest atoms to subselection, not including the 
   * subselection itself.  If there are only j<k atoms not in subselection,
   * then only j atoms will be selected. */
  PredicatePtr k_nearest_predicate(SystemPtr sys, unsigned k, PredicatePtr S);

}}}

#endif
