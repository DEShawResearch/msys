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
    table->category = schema->category;
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
    std::set<std::string> TableSchemas() {
        std::set<std::string> s;
        for (unsigned i=0; i<schema_count(); i++) {
            s.insert(schema_name(i));
        }
        return s;
    }
    std::set<std::string> NonbondedSchemas() {
        std::set<std::string> s;
        for (unsigned i=0; i<nonbonded_schema_count(); i++) {
            s.insert(nonbonded_schema_name(i));
        }
        return s;
    }

    TermTablePtr AddTable(SystemPtr sys, const std::string& type,
                                         const std::string& _name) {

        std::string name = _name.size() ? _name : type;
        TermTablePtr table = sys->table(name);
        if (table) return table;
        const schema_t* schema = find_schema(type);
        if (!schema) {
            std::stringstream ss;
            ss << "Unknown DMS table schema '" << type << "'";
            throw std::runtime_error(ss.str());
        }
        table = sys->addTable(name, schema->nsites);
        configure_table(schema, table);
        return table;
    }

    TermTablePtr AddNonbonded(SystemPtr sys, 
                                 const std::string& funct,
                                 const std::string& rule) {
        const schema_t* schema = find_nonbonded_schema(funct);
        if (!schema) {
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
        configure_table(schema, table);
        return table;
    }
    
}}
