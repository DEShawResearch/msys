#ifndef msys_schema_hxx
#define msys_schema_hxx

#include "system.hxx"

namespace desres { namespace msys {

    /* Names of available table schemas. */
    std::set<std::string> TableSchemas();

    /* Names of available nonbonded schemas. */
    std::set<std::string> NonbondedSchemas();

    /* Add a table to the system if it not already present,
     * returning it.  If optional name field is provided, the table
     * will be added with the given name; otherwise the name is taken
     * from the table schema. */
    TermTablePtr AddTable(SystemPtr sys, const std::string& type,
                                         const std::string& name="");

    /* Add a nonbonded table to the system, and configure the nonbonded
     * info according to funct and rule.  funct must be the name of recognized
     * nonbonded type.  rule is not checked; at some point in the future we
     * might start requiring that it be one of the valid combining rules for
     * the specified funct.  If nonbonded_info's vdw_funct and vdw_rule 
     * are empty, they are overridden by the provided values; otherwise, the
     * corresponding values must agree if funct and rule are not empty.
     * A nonbonded table is returned.  */
    TermTablePtr AddNonbonded(SystemPtr sys, 
                              const std::string& funct,
                              const std::string& rule);
}}

#endif
