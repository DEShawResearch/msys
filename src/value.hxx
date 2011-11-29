#ifndef mol_value_hxx
#define mol_value_hxx

#include "types.hxx"

namespace desres { namespace msys {

    enum ValueType {
        IntType,
        FloatType,
        StringType
    };
    
    union Value {
        Int     i;
        Float   f;
        char *  s;
    };
    
    class ValueRef {
        ValueType   _type;
        Value&      _val;
    
    public:
        ValueRef( ValueType type, Value& val ) : _type(type), _val(val) {}
        ValueType type() const { return _type; }
    
        Int asInt() const;
        Float asFloat() const;
        String asString() const;
    
        void fromInt(const Int& i);
        void fromFloat(const Float& i);
        void fromString(const String& i);
    
        void assign(uint32_t v){ fromInt(v); }
        void assign(short v)   { fromInt(v); }
        void assign(int32_t v) { fromInt(v); }
        void assign(int64_t v) { fromInt(v); }
    
        void assign(float v)  { fromFloat(v); }
        void assign(double v) { fromFloat(v); }
    
        void assign(const std::string& v) { fromString(v); }
        /* don't die on NULL - treat as empty string */
        void assign(const char * s ) { fromString( s ? s : "" ); }

        void assign(const ValueRef& v) {
            switch (type()) {
                default:
                case IntType:    assign(v.asInt()); break;
                case FloatType:  assign(v.asFloat()); break;
                case StringType: assign(v.asString()); break;
            }
        }
    
        /* int conversions allowed from IntType values */
        operator short() const { return asInt(); }
        operator int() const   { return asInt(); }
        operator long() const  { return asInt(); }
    
        /* float conversions allowed from FloatType values */
        operator float() const  { return asFloat(); }
        operator double() const { return asFloat(); }
    
        /* string conversions allowed from StringType values */
        operator std::string() const { return asString(); }
    
        /* comparision to a concrete type is performed by converting this
         * to the type of the argument. */
        template <typename T>
        bool operator==(const T& rhs) const {
            T lhs = *this;
            return lhs==rhs;
        }

        /* comparison between IntType and FloatType invokes the 
         * underlying C comparison of Int and Float.  Comparison of strings
         * uses strcmp.  Any other comparison returns false. */
        bool operator==(const ValueRef& rhs) const;

        /* return -1, 0 or 1 if this compares less than, equal to, or
         * greater than rhs */
        int compare(const ValueRef& rhs) const;

        template <typename T>
        bool operator!=(const T& rhs) const {
            return ! operator==(rhs);
        }
    
        /* assignment */
        template <typename T>
        ValueRef& operator=(const T& rhs) {
            assign(rhs);
            return *this;
        }
    
        ValueRef& operator=(const ValueRef& rhs) {
            assign(rhs);
            return *this;
        }
    };

}}

#endif
