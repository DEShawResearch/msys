#ifndef msys_atomsel_within_predicate_hxx
#define msys_atomsel_within_predicate_hxx

#include "predicate.hxx"
#include "../system.hxx"

namespace desres { namespace msys { namespace atomsel {

  /*! true for elements whose cartesian distance from elements in the 
   * subselection is less than r */
  PredicatePtr within_predicate( SystemPtr sys, double r, PredicatePtr S );

  /*! true for elements whose cartesian distance from elements in the 
   * subselection is less than r, and are also not in S */
  PredicatePtr exwithin_predicate( SystemPtr sys, double r, PredicatePtr S );

  /*! true for elements whose minimum image distance from elements in the 
   * subselection is less than r */
  PredicatePtr pbwithin_predicate( SystemPtr sys, double r, PredicatePtr S );

  /*! true for elements within n bonds of subselection */
  PredicatePtr withinbonds_predicate( SystemPtr sys, int n, PredicatePtr S );

  /*! true for the k nearest atoms to subselection, not including the 
   * subselection itself.  If there are only j<k atoms not in subselection,
   * then only j atoms will be selected. */
  PredicatePtr k_nearest_predicate(SystemPtr sys, unsigned k, PredicatePtr S);
  PredicatePtr k_pbnearest_predicate(SystemPtr sys, unsigned k, PredicatePtr S);

}}}

#endif
