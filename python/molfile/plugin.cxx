#include "molfilemodule.hxx"
#include <vector>
#include <map>
#include <algorithm>

using namespace desres::molfile;

namespace {

    bool compare_bond( const bond_t& b1, const bond_t& b2) {
        if (b1.from!=b2.from) return b1.from<b2.from;
        return b1.to<b2.to;
    }

    /* convert from Python representation of bonds/atoms to the C++ rep */
    int assemble_atoms( list atoms, std::vector<atom_t>& atomlist,
                                    std::vector<bond_t>& bondlist ) {
        int optflags=0;
        if (atoms.is_none()) {
            atomlist.resize(0);
            bondlist.resize(0);
            return 0;
        }

        /* first pass: hash the atom objects */
        typedef std::map<PyObject *, Py_ssize_t> AtomHash;
        AtomHash atom_hash;
        atomlist.resize(len(atoms));
        for (size_t i=0; i<len(atoms); i++) {
            handle obj = atoms[i];
            PyObject * ptr = obj.ptr();
            if (ptr->ob_type != &AtomType) {
                PyErr_Format(PyExc_ValueError,
                        "atoms must contain only Atom items");
                throw error_already_set();
            }
            Atom_t *atom = reinterpret_cast<Atom_t*>(ptr);
            atomlist[i] = atom->atom;
            optflags |= atom->optflags;
            atom_hash[ptr]=i;
        }

        /* second pass: assemble the bonds */
        for (size_t i=0; i<len(atoms); i++) {
            object obj=atoms[i];
            Atom_t *atom = reinterpret_cast<Atom_t*>(obj.ptr());
            /* iterate over keyvals in bond dict */
            Py_ssize_t pos=0;
            PyObject * key, * val;
            while (PyDict_Next(atom->bonds, &pos, &key, &val)) {
                AtomHash::const_iterator j=atom_hash.find(key);
                if (j!=atom_hash.end() && ssize_t(i)<j->second) {
                    float order = PyFloat_AsDouble(val);
                    if (PyErr_Occurred()) throw error_already_set();
                    bondlist.push_back(bond_t(i,j->second,order));
                }
            }
        }
        std::sort(bondlist.begin(), bondlist.end(), compare_bond);
        return optflags;
    }

    Writer * plugin_write( const molfile_plugin_t& self, 
                           const std::string& path,
                           object& atomsobj,
                           object& natomsobj ) {


        /* able to write? */
        if (!self.open_file_write) {
            PyErr_Format(PyExc_ValueError, "Plugin does not support writing");
            throw error_already_set();
        }

        /* get natoms from either atoms list or specified natoms */
        ssize_t natoms=-1;
        if (!atomsobj.is_none())       natoms = len(atomsobj);
        else if (!natomsobj.is_none()) natoms = natomsobj.cast<ssize_t>();
        else if (self.write_structure) {
            PyErr_Format(PyExc_ValueError, "Provide either atoms or natoms>0");
            throw error_already_set();
        }

        Writer * writer = new Writer(&self, path.c_str(), natoms);
        if (self.write_structure) try {
            std::vector<atom_t> atoms;
            std::vector<bond_t> bonds;
            int optflags = assemble_atoms( atomsobj, atoms, bonds );
            writer->write_atoms( atoms, bonds, optflags );
        } catch (std::exception& e) {
            delete writer;
            PyErr_Format(PyExc_RuntimeError, "Writing structure failed: %s",
                    e.what());
            throw error_already_set();
        }

        return writer;
    }

    tuple plugin_version(const molfile_plugin_t& p) {
        return make_tuple(p.majorv, p.minorv);
    }
    str plugin_extensions(const molfile_plugin_t& p) { 
        const char * ext = p.filename_extension;
        return ext ? ext : "";
    }
    bool plugin_can_read(const molfile_plugin_t& p) {
        return p.open_file_read != NULL;
    }
    bool plugin_can_write(const molfile_plugin_t& p) {
        return p.open_file_write != NULL;
    }

    object plugin_repr(const molfile_plugin_t& p) {
        return str("<Plugin for %s>").format(p.prettyname);
    }

    Reader * plugin_read(const molfile_plugin_t& p, 
                         const std::string& path,
                         bool double_precision) {
        return new Reader(&p, path.c_str(), double_precision);
    }
}

void desres::molfile::export_plugin(module m) {

    class_<molfile_plugin_t>(m, "Plugin", "Interface to readers and writers")
        .def(init<>())
        .def_readonly("name", &molfile_plugin_t::name)
        .def_readonly("prettyname", &molfile_plugin_t::prettyname)
        .def_property_readonly("version", plugin_version)
        .def_property_readonly("filename_extensions", plugin_extensions)
        .def_property_readonly("can_read", plugin_can_read)
        .def_property_readonly("can_write", plugin_can_write)
        .def("__repr__", plugin_repr)
        .def("read", plugin_read, 
                return_value_policy::reference,
                 arg("path")
                ,arg("double_precision")=false,
                "Open a file for reading")
        .def("write", plugin_write,
                return_value_policy::reference,
                 arg("path"),
                 arg("atoms")=none(),
                 arg("natoms")=none(),
                 "write(path,atoms=None,natoms=None)")
        ;
}

