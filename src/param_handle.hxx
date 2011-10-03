#ifndef msys_param_handle_hxx
#define msys_param_handle_hxx

#include "param_table.hxx"
#include "ref.hxx"
#include "append.hxx"

namespace desres { namespace msys {

    class ParamTableHandle;

    /* a Param refers to a row in the table of parameters */
    class Param {
        ParamTable* _params;
        Id      _id;
    
    public:
        Param( ParamTable* params=NULL, Id id=BadId ) 
        : _params(params), _id(id) {}

        inline ParamTableHandle params();
        
        bool operator==(const Param& other) const {
            return _params==other._params && _id==other._id;
        }

        Id id() const { return _id; }
        ValueType type(Id i) const { return _params->propType(i); }
        ValueType type(const String& name) const { 
            return type(_params->propIndex(name));
        }
        String name(Id i) const { return _params->propName(i); }
        ValueRef val(Id i)      { return _params->value(id(), i); }
        ValueRef val(const String& name) {
            return val(_params->propIndex(name));
        }
    };

    class ParamTableHandle {
        ParamTablePtr _ptr;

    public:
        /* default constructor creates invalid handle */
        ParamTableHandle() {}

        /* convertible from ParamTablePtr */
        ParamTableHandle(ParamTablePtr ptr) : _ptr(ptr) {}

        ParamTablePtr ptr() const { return _ptr; }

        bool operator==(const ParamTableHandle& other) const {
            return _ptr==other._ptr;
        }
        bool operator!=(const ParamTableHandle& other) const {
            return _ptr!=other._ptr;
        }

        Id paramCount() const { return _ptr->paramCount(); }
        Param param(Id id) const { return Param(_ptr.get(), id); }
        Param addParam() {
            return Param(_ptr.get(), _ptr->addParam());
        }

        Id propCount() const { return _ptr->propCount(); }
        String propName(Id i) const { return _ptr->propName(i); }
        ValueType propType(Id i) const { return _ptr->propType(i); }
        Id propIndex(const String& name) const { return _ptr->propIndex(name); }
        Id addProp(const String& name, ValueType type) {
            return _ptr->addProp(name,type);
        }

        IdList appendParams( ParamTableHandle src, 
                             IdList const& params ) {
            return AppendParams( ptr(), src.ptr(), params );
        }


    };

    ParamTableHandle Param::params() {
        return _params->shared_from_this();
    }

}}

#endif
