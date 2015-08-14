#define NO_IMPORT_ARRAY 1
#include "unique_symbol.hxx"
#include <boost/python.hpp>
#include "molecule.hxx"

using namespace boost::python;

namespace {
    using namespace desres::msys;
    Molecule* molecule_iterator_next(MoleculeIterator& self) {
        return self.next().release();
    }
    str mol_getitem(Molecule const& m, std::string const& key) {
        return str(m.data().at(key));
    }
    void mol_setitem(Molecule& m, std::string const& key, std::string val) {
        m.data()[key] = val;
    }
    void mol_delitem(Molecule& m, std::string const& key) {
        m.data().erase(key);
    }

    object mol_get(Molecule const& m, std::string const& key, object defval) {
        auto it = m.data().find(key);
        if (it==m.data().end()) return defval;
        return str(it->second);
    }

    std::string get_name(Molecule const& m) {
        return m.name();
    }
    void set_name(Molecule& m, std::string const& s) {
        m.name() = s;
    }

    list mol_keys(Molecule const& m) {
        list L;
        for (auto const& it : m.data()) {
            L.append(object(it.first));
        }
        return L;
    }
}

namespace desres { namespace msys {

    void export_molecule() {
        class_<Molecule, boost::noncopyable>("Molecule", no_init)
            .add_property("name",   get_name, set_name)
            .add_property("natoms", &Molecule::natoms)
            .add_property("nbonds", &Molecule::nbonds)
            .def("__getitem__", mol_getitem)
            .def("__setitem__", mol_setitem)
            .def("__delitem__", mol_delitem)
            .def("get",         mol_get,
                    (arg("key"),
                     arg("defval")=object()))
            .def("keys",        mol_keys)
        ;

        class_<MoleculeIterator, boost::noncopyable>("MoleculeIterator", no_init)
            .def("next", molecule_iterator_next,
                    return_value_policy<manage_new_object>())
            ;
    }
}}

