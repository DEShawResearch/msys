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

        /* number of TermTables holding this ParamTable */
        Id _refcnt;
    
        ParamTable();
    public:
        static boost::shared_ptr<ParamTable> create();
        ~ParamTable();

        /* Do multiple TermTables use this ParamTable? */
        bool shared() const { return _refcnt > 1; }

        /* incref and decref should be used only from within TermTable
         * to indicate whether a ParamTable is shared by multiple 
         * term tables. */
        void incref() { ++_refcnt; }
        void decref() { --_refcnt; }
    
        Id paramCount() const { return _nrows; }
        Id addParam();
        bool hasParam(Id param) const {
            return param<paramCount();
        }
        /* create a duplicate of the given parameter set, returning
         * its id.  If param is BadId, default parameters are used;
         * this is identical to addParam(). */
        Id duplicate(Id param);

    
        Id propCount() const            { return _props.size();     }
        String propName(Id i) const     { return _props.at(i).name; }
        ValueType propType(Id i) const  { return _props.at(i).type; }
        Id propIndex(const String& name) const;
        Id addProp(const String& name, ValueType type);

        /* delete the property with the given index.  If index is BadId,
         * no action is taken; otherwise it is an error if the index is
         * invalid. */
        void delProp(Id index);

        ValueRef value(Id row, Id col)  { 
            return ValueRef(propType(col), _props.at(col).vals.at(row)); 
        }
        ValueRef value(Id row, String const& name)  { 
            Id col = propIndex(name);
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
