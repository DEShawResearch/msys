#include "wrap_obj.hxx"
#include "value.hxx"
#include <boost/foreach.hpp>

using namespace boost::python;
using namespace desres::msys;


namespace {

    object propmap_keys(VariantMap& p) {
        list s;
        BOOST_FOREACH(VariantMap::value_type keyvals, p) {
            s.append(keyvals.first);
        }
        return s;
    }

    class get_visitor : public boost::static_visitor<> {
        object& o;
    public:
        explicit get_visitor(object& o) : o(o) {}
        void operator()(Int& i) const { o=object(i); }
        void operator()(Float& i) const { o=object(i); }
        void operator()(String& i) const { o=object(i); }
    };

    object propmap_get(VariantMap& p, std::string const& key) {
        VariantMap::iterator i=p.find(key);
        if (i==p.end()) return object();
        object o;
        boost::apply_visitor(get_visitor(o), i->second);
        return o;
    }

    void propmap_set(VariantMap& p, std::string const& key, object type, 
                     object val) {
        ValueType t = as_value_type(type);
        Variant& v = p[key];
        switch (t) {
            case IntType: v = extract<int64_t>(val); break;
            case FloatType: v = extract<double>(val); break;
            case StringType:
            default:        v = extract<std::string>(val); break;
        }
    }
    void propmap_del(VariantMap& p, std::string const& key) {
        VariantMap::iterator i=p.find(key);
        if (i!=p.end()) p.erase(i);
    }
}

namespace desres { namespace msys { 

    void export_propmap() {
        class_<VariantMap>("VariantMap", no_init)
            .def("empty", &VariantMap::empty)
            .def("keys",  propmap_keys)
            .def("get",   propmap_get)
            .def("set",   propmap_set)
            .def("_del",  propmap_del)
            ;
    }

}}

