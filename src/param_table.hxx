#ifndef msys_param_table_hxx
#define msys_param_table_hxx

#include "value.hxx"
#include <deque>
#include <map>
#include <boost/enable_shared_from_this.hpp>

namespace desres { namespace msys {

    /* A ParamTable instance is a table whose rows are indexed by id and whose
     * columns are named and can have int, float, and string types. */
    class ParamTable : public boost::enable_shared_from_this<ParamTable> {
    
        /* we use deque so that references to values returned by the value()
         * method aren't invalidated by later calls to addProp(). */
        typedef std::deque<Value> ValueList;
    
        struct Property {
            String      name;
            ValueType   type;
            ValueList   vals;
    
            void extend();
        };
    
        typedef std::deque<Property> PropList;
        PropList _props;
    
        /* number of rows in each property */
        Id  _nrows;
    
        ParamTable();
    public:
        static boost::shared_ptr<ParamTable> create();
        ~ParamTable();
    
        Id paramCount() const { return _nrows; }
        Id addParam();
        bool hasParam(Id param) const {
            return param<paramCount();
        }
    
        Id propCount() const            { return _props.size();     }
        String propName(Id i) const     { return _props.at(i).name; }
        ValueType propType(Id i) const  { return _props.at(i).type; }
        Id propIndex(const String& name) const;
        Id addProp(const String& name, ValueType type);
        ValueRef value(Id row, Id col)  { 
            return ValueRef(propType(col), _props.at(col).vals.at(row)); 
        }
        IdList params() const {
            IdList p(_nrows);
            for (Id i=0; i<_nrows; i++) p[i]=i;
            return p;
        }
    };
    typedef boost::shared_ptr<ParamTable> ParamTablePtr;
}}

#endif
