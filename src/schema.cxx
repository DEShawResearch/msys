#include "schema.hxx"
#include "schema/schema.hxx"
#include "term_table.hxx"

#include <sstream>
#include <stdio.h>
#include <string.h>
#include <stdexcept>

#include <boost/algorithm/string.hpp> /* for boost::trim */

using namespace desres::msys;

static void configure_table(const schema_t* schema, 
                            TermTablePtr table) {
    table->category = parse(schema->category);
    const column_t* props;
    props = schema->param_props;
    for (int i=0; i<DMS_MAX_PARAM_PROPS; i++) {
        if (!props[i].name) break;
        ValueType type = ValueType(props[i].type);
        table->params()->addProp(props[i].name, type);
    }
    props = schema->term_props;
    for (int i=0; i<DMS_MAX_TERM_PROPS; i++) {
        if (!props[i].name) break;
        ValueType type = ValueType(props[i].type);
        table->addTermProp(props[i].name, type);
    }
}

namespace desres { namespace msys {
    std::vector<std::string> TableSchemas() {
        std::vector<std::string> s;
        for (unsigned i=0; i<schema_count(); i++) {
            s.push_back(schema(i).name);
        }
        return s;
    }
    std::vector<std::string> NonbondedSchemas() {
        std::vector<std::string> s;
        for (unsigned i=0; i<nonbonded_schema_count(); i++) {
            s.push_back(nonbonded_schema(i).name);
        }
        return s;
    }

    TermTablePtr AddTable(SystemPtr sys, const std::string& type,
                                         const std::string& _name) {

        std::string name = _name.size() ? _name : type;
        TermTablePtr table = sys->table(name);
        if (table) return table;
        for (unsigned i=0; i<schema_count(); i++) {
            schema_t s = schema(i);
            if (type!=s.name) continue;
            table = sys->addTable(name, s.nsites);
            configure_table(&s, table);
        }
        if (!table) MSYS_FAIL("Unknown DMS table schema: " << type); 
        return table;
    }

    TermTablePtr AddNonbonded(SystemPtr sys, 
                                 const std::string& funct,
                                 const std::string& rule) {
        schema_t s;
        for (unsigned i=0; i<nonbonded_schema_count(); i++) {
            s=nonbonded_schema(i);
            if (funct==s.name) break;
            s.name=NULL;
        }
        if (!s.name) {
            std::stringstream ss;
            ss << "Unknown DMS nonbonded table type '" << funct << "'";
            throw std::runtime_error(ss.str());
        }
        NonbondedInfo& info = sys->nonbonded_info;
        if (rule.size()) {
            std::string newrule(rule), oldrule(info.vdw_rule);
            boost::to_lower(newrule);
            boost::to_lower(oldrule);
            if (oldrule.size() && newrule != oldrule) {
                std::stringstream ss;
                ss << "Cannot add nonbonded table with vdw_rule '" << rule
                    << "' because nonbonded_info.vdw_rule is already set to '"
                    << info.vdw_rule << "'";
                throw std::runtime_error(ss.str());
            }
            info.vdw_rule = newrule;
        }
        if (funct.size()) {
            std::string newfunct(funct), oldfunct(info.vdw_funct);
            boost::to_lower(newfunct);
            boost::to_lower(oldfunct);
            if (oldfunct.size() && newfunct!= oldfunct) {
                std::stringstream ss;
                ss << "Cannot add nonbonded table with vdw_funct'" << funct
                    << "' because nonbonded_info.vdw_funct is already set to '"
                    << info.vdw_funct;
                throw std::runtime_error(ss.str());
            }
            info.vdw_funct = newfunct;
        }
        TermTablePtr table = sys->table("nonbonded");
        if (table) return table;
        table = sys->addTable("nonbonded", 1);
        configure_table(&s, table);
        return table;
    }
    
}}
