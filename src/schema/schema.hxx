#ifndef msys_dms_dms_schema_hxx
#define msys_dms_dms_schema_hxx

#include "../value.hxx"

#define DMS_MAX_TERM_PROPS  6
#define DMS_MAX_PARAM_PROPS 20

namespace desres { namespace msys { 

    struct column_t {
        const char* name;
        int         type;
        const char* comment;
        column_t(const char *_name=NULL, int _type=1, const char *_comment="") 
        : name(_name), type(_type), comment(_comment) {}
    };


    struct schema_t {
        const char* name;
        const char* category;
        int         nsites;
        column_t param_props[DMS_MAX_PARAM_PROPS];
        column_t term_props[DMS_MAX_TERM_PROPS];
    };

    unsigned schema_count();
    const char* schema_name(unsigned i);

    unsigned nonbonded_schema_count();
    const char* nonbonded_schema_name(unsigned i);

    const schema_t* find_schema(const std::string& name);
    const schema_t* find_nonbonded_schema(const std::string& name);


}}

#endif
