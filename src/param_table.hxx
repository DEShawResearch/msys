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
    
        struct Property : public ValueCallback {
            String      name;
            ValueType   type;
            ValueList   vals;
    
            /* index on values in this column */
            typedef std::map<ValueRef, IdList> Index;
            Index       index;  

            /* one more than last row indexed */
            Id          maxIndexId;
            
            /* reindex rows [maxIndexId, vals.size()) */
            void update_index();

            /* add a new row */
            void extend();

            /* invalidate the entire index */
            virtual void valueChanged() {
                index.clear();
                maxIndexId = 0;
            }

            Property() : maxIndexId(0) {}
        };
        
        typedef std::deque<Property> PropList;
        PropList _props;
    
        /* number of rows in each property */
        Id  _nrows;

        /* reference count for parameters in this table */
        IdList _paramrefs;

        ParamTable();
    public:
        static boost::shared_ptr<ParamTable> create();
        ~ParamTable();

        /* increment/decrement the reference count of the given parameter */
        void incref(Id param);
        void decref(Id param);

        /* reference for the given parameter */
        Id refcount(Id param) const {
            return param>=_paramrefs.size() ? 0 : _paramrefs[param];
        }
    
        Id paramCount() const { return _nrows; }
        Id addParam();
        bool hasParam(Id param) const {
            return param<paramCount();
        }
        /* create a duplicate of the given parameter set, returning
         * its id.  If param is BadId, default parameters are used;
         * this is identical to addParam(). */
        Id duplicate(Id param);

        /* Return -1, 0, 1 if the values of L are less than, equal to,
         * or greater than those of R, with comparisons performed
         * lexicographically using ValueRef::compare(). */
        int compare(Id L, Id R);
    
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
            return ValueRef(
                    propType(col), 
                    _props.at(col).vals.at(row),
                    &_props.at(col));
        }
        ValueRef value(Id row, String const& name);

        IdList params() const {
            IdList p(_nrows);
            for (Id i=0; i<_nrows; i++) p[i]=i;
            return p;
        }

        /* find parameters with the given value in the given column */
        IdList findInt(Id col, Int const& val);
        IdList findFloat(Id col, Float const& val);
        IdList findString(Id col, String const& val);
        IdList findValue(Id col, ValueRef const& v);
    };
    typedef boost::shared_ptr<ParamTable> ParamTablePtr;
}}

#endif
