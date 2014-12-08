#ifndef desres_msys_prop_map_hxx
#define desres_msys_prop_map_hxx

#include "value.hxx"
#include <map>

namespace desres { namespace msys {

    class PropertyMap {

        typedef std::pair<ValueType, Value> Property;
        typedef std::map<String, Property> Map;

        Map _map;

    public:
        ~PropertyMap();

        /* anything here? */
        bool empty() const { return _map.empty(); }

        /* list of keys */
        std::vector<String> keys() const;

        /* type for given key */
        ValueType type(String const& key) const;

        /* is key present? */
        bool has(String const& key) const;

        /* add a property with a given type */
        ValueRef add(String const& key, ValueType type);

        /* remove a property.  Exception if missing */
        void del(String const& key);

        /* get value for read/write */
        ValueRef get(String const& key);
    };
}}

#endif

