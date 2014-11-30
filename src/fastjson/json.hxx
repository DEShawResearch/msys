/* @COPYRIGHT@ */

#ifndef desres_fastjson_hxx
#define desres_fastjson_hxx

#include <cstring>

namespace desres { namespace msys { namespace fastjson {

    struct keyval_t;

    class Json {

    public:
        enum kind_t {
            Invalid,
            Array,
            Object,
            Int,
            Float,
            String,
            Null,
            Bool
        };

    private:
        struct array_t {
            Json * elems;
            int num;
            int max;
        };

        struct object_t {
            keyval_t * elems;
            int num;
            int max;
        };

        kind_t _kind;

        union {
            int      i;
            double   f;
            char *   s;
            array_t  a;
            object_t o;
        } u;

        /* no copy */
        Json(const Json&) {}
        Json& operator=(const Json&) { return *this; }

        /* release internal memory */
        void clear();

    public:

        /* default constructor has no type */
        Json() : _kind(Invalid) {}

        /* destructor */
        ~Json() { clear(); }

        /* become an empty array, releasing any previous data */
        Json& to_array();

        /* become an empty object, releasing any previous data */
        Json& to_object();

        /* become a json of null type, releasing any previous data */
        Json& to_null() {
            clear();
            _kind = Null;
            return *this;
        }

        /* become other primitive types */
        Json& to_int(int v) {
            clear();
            u.i=v;
            _kind = Int;
            return *this;
        }

        Json& to_float(double v) {
            clear();
            u.f=v;
            _kind = Float;
            return *this;
        }

        Json& to_string(const char * s) {
            clear();
            u.s=strdup(s);
            _kind = String;
            return *this;
        }

        Json& to_bool(bool v) {
            clear();
            u.i=v;
            _kind = Bool;
            return *this;
        }

        /* get int, float, string, or bool representations.
         *
         * int is castable to double and bool
         * bool is castable to int 
         * Any other conversion attempt will throw std::invalid_argument.
         */

        int          as_int() const;
        double       as_float() const;
        const char * as_string() const;
        bool         as_bool() const;

        /* get representations, supplying a default value if json is invalid 
         * or null */
        int as_int(int defval) const {
            if (_kind == Null || _kind == Invalid ) return defval;
            return as_int();
        }
        double as_float(double defval) const {
            if (_kind == Null || _kind == Invalid ) return defval;
            return as_float();
        }
        const char * as_string(const char * defval) const {
            if (_kind == Null || _kind == Invalid ) return defval;
            return as_string();
        }
        bool as_bool(bool defval) const {
            if (_kind == Null || _kind == Invalid ) return defval;
            return as_bool();
        }

        /* what kind of json? */
        kind_t kind() const { return _kind; }

        /* string representation of the kind */
        const char * kindstr() const;

        /* is json valid/invalid? */
        bool operator!() const { return _kind == Invalid; }
        bool valid()     const { return _kind != Invalid; }

        /* number of elements, if array or object. 
         * If not one of these kinds, throw std::invalid_argument.
         */
        int size() const;

        /* nth array or object element, for modification */
        Json& elem(int n);

        /* nth array or object element.  If this is invalid json, returns
         * invalid json. */
        const Json& elem(int n) const;

        /* nth key, if object */
        const char * key(int n) const;

        /* value for key, if object.  Returns invalid json on nonexistent key,
         * or when called on invalid json */
        const Json& get(const char * key) const;

        /* swap storage from the given json. O(1) complexity.  */
        void swap(Json& js);

        /* copy data from the given json.  O(sizeof(js)) complexity. */
        Json& copy(const Json& js);

        /* recursively test for equality.  
         *
         * Jsons of Object type are compared without respect to key
         * order; i.e., jsons with the same keys in a different order
         * will compare equal.  */
        bool equals(const Json& js) const;

        /* append json to array, stealing its storage. */
        Json& append( Json& js );

        /* add keyval to array, stealing storage of val */
        Json& append( const char * key, Json& val );

    };

    struct keyval_t {
        char * key;
        Json   val;
    };

}}}

#endif

