#include "molfilemodule.hxx"
#include <boost/python.hpp>
#include <vector>
#include <map>


using namespace desres::molfile;
using namespace boost::python;

namespace {

    bool compare_bond( const bond_t& b1, const bond_t& b2) {
        if (b1.from!=b2.from) return b1.from<b2.from;
        return b1.to<b2.to;
    }

    /* convert from Python representation of bonds/atoms to the C++ rep */
    int assemble_atoms( object& atoms, std::vector<atom_t>& atomlist,
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
        for (Py_ssize_t i=0; i<len(atoms); i++) {
            object obj = atoms[i];
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
        for (Py_ssize_t i=0; i<len(atoms); i++) {
            object obj=atoms[i];
            Atom_t *atom = reinterpret_cast<Atom_t*>(obj.ptr());
            /* iterate over keyvals in bond dict */
            Py_ssize_t pos=0;
            PyObject * key, * val;
            while (PyDict_Next(atom->bonds, &pos, &key, &val)) {
                AtomHash::const_iterator j=atom_hash.find(key);
                if (j!=atom_hash.end() && i<j->second) {
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
        Py_ssize_t natoms=-1;
        if (!atomsobj.is_none())       natoms = len(atomsobj);
        else if (!natomsobj.is_none()) natoms = extract<Py_ssize_t>(natomsobj);
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
        return "<Plugin for %s>" % make_tuple(p.prettyname);
    }

    Reader * plugin_read(const molfile_plugin_t& p, 
                         const std::string& path,
                         bool double_precision,
                         bool with_gids) {
        return new Reader(&p, path.c_str(), double_precision, with_gids);
    }
}

void desres::molfile::export_plugin() {

    class_<molfile_plugin_t>("Plugin", init<>())
        .add_property("name", &molfile_plugin_t::name)
        .add_property("prettyname", &molfile_plugin_t::prettyname)
        .add_property("version", plugin_version)
        .add_property("filename_extensions", plugin_extensions)
        .add_property("can_read", plugin_can_read)
        .add_property("can_write", plugin_can_write)
        .def("__repr__", plugin_repr)
        .def("read", plugin_read, 
                (arg("path")
                ,arg("double_precision")=false
                ,arg("with_gids")=false),
                return_value_policy<manage_new_object,
                return_internal_reference<1> >())
        .def("write", plugin_write,
                (arg("path"),
                 arg("atoms")=object(),
                 arg("natoms")=object()),
                 "write(path,atoms=None,natoms=None)",
                return_value_policy<manage_new_object,
                return_internal_reference<1> >())
        ;
}

