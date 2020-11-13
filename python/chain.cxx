#include "pymod.hxx"
#include <pybind11/stl.h>
#include <pybind11/operators.h>
#include <msys/system.hxx>

namespace desres { namespace msys { 

    template <typename T>
    struct Handle {
        SystemPtr mol;
        Id id;

        bool operator==(T const& rhs) const { return mol==rhs.mol && id==rhs.id; }
        bool operator!=(T const& rhs) const { return mol!=rhs.mol || id!=rhs.id; }
    };

    struct Atom : Handle<Atom> { Atom(SystemPtr m, Id i) : Handle{m,i} {} };
    struct Bond : Handle<Bond> { Bond(SystemPtr m, Id i) : Handle{m,i} {} };
    struct Residue : Handle<Residue> {
        Residue(SystemPtr m, Id i) : Handle{m,i} {}
        residue_t& data() { return mol->residue(id); }
        const char* name() { return data().name.c_str(); }
        void setName(std::string const& s) { data().name=s; }
        const char* insertion() { return data().insertion.c_str(); }
        void setInsertion(std::string const& s) { data().insertion=s; }
        int resid() { return data().resid; }
        void setResid(int i) { data().resid=i; }
        void remove() { mol->delResidue(id); }
        Id addAtom() { return mol->addAtom(id); }
        IdList atoms() { return mol->atomsForResidue(id); }
        Id natoms() { return mol->atomCountForResidue(id); }
        Id chainId() { return data().chain; }
        void setChainId(Id chn) { mol->setChain(id, chn); }
    };

    struct Chain : Handle<Chain> {
        Chain(SystemPtr m, Id i) : Handle{m,i} {}
        chain_t& data() { return mol->chain(id); }
        const char* name() { return data().name.c_str(); }
        void setName(std::string const& s) { data().name=s; }
        const char* segid() { return data().segid.c_str(); }
        void setSegid(std::string const& s) { data().segid=s; }
        Id ctId() { return data().ct; }
        void setCtId(Id ct) { mol->setCt(id, ct); }
        IdList residues() { return mol->residuesForChain(id); }
        Id nresidues() { return mol->residueCountForChain(id); }
        Id addResidue() { return mol->addResidue(id); }
        void remove() { mol->delChain(id); }
    };

    void export_chain(module m) {

        class_<Atom>(m, "Atom")
        .def(init<SystemPtr, Id>())
        .def(self == self)
        .def(self != self)
        .def("__hash__", [](Atom const& a) { return hash(make_tuple(a.mol.get(), a.id)); })
        .def("__lt__", [](Atom& lhs, Atom& rhs) { return (lhs.mol < rhs.mol) || (lhs.id < rhs.id); })
        .def_readonly("_ptr", &Atom::mol)
        .def_readonly("id", &Atom::id)
        .def_readonly("_id", &Atom::id) // backward compatability
        .def_property("x", [](Atom& a) { return a.mol->atom(a.id).x; },
                [](Atom& a, double val) { a.mol->atom(a.id).x = val; },
                "x component of position")
        .def_property("y", [](Atom& a) { return a.mol->atom(a.id).y; },
                [](Atom& a, double val) { a.mol->atom(a.id).y = val; },
                "y component of position")
        .def_property("z", [](Atom& a) { return a.mol->atom(a.id).z; },
                [](Atom& a, double val) { a.mol->atom(a.id).z = val; },
                "z component of position")
        .def_property("vx", [](Atom& a) { return a.mol->atom(a.id).vx; },
                [](Atom& a, double val) { a.mol->atom(a.id).vx = val; },
                "x component of velocity")
        .def_property("vy", [](Atom& a) { return a.mol->atom(a.id).vy; },
                [](Atom& a, double val) { a.mol->atom(a.id).vy = val; },
                "y component of velocity")
        .def_property("vz", [](Atom& a) { return a.mol->atom(a.id).vz; },
                [](Atom& a, double val) { a.mol->atom(a.id).vz = val; },
                "z component of velocity")
        .def_property("mass", [](Atom& a) { return a.mol->atom(a.id).mass; },
                [](Atom& a, double val) { a.mol->atom(a.id).mass = val; },
                "mass")
        .def_property("charge", [](Atom& a) { return a.mol->atom(a.id).charge; },
                [](Atom& a, double val) { a.mol->atom(a.id).charge= val; },
                "charge")
        .def_property("atomic_number", [](Atom& a) { return a.mol->atom(a.id).atomic_number; },
                [](Atom& a, int val) { a.mol->atom(a.id).atomic_number = val; },
                "atomic number")
        .def_property("formal_charge", [](Atom& a) { return a.mol->atom(a.id).formal_charge; },
                [](Atom& a, int val) { a.mol->atom(a.id).formal_charge = val; },
                "formal charge")
        .def_property("name", [](Atom& a) { return str(a.mol->atom(a.id).name); },
                [](Atom& a, std::string const& val) { a.mol->atom(a.id).name = val; },
                "name")
        .def_property_readonly("fragid", [](Atom& a) { return a.mol->atom(a.id).fragid; },
                "fragment id")
        .def_property_readonly("resid", [](Atom& a) { return a.mol->atom(a.id).residue; },
                "residue id")
        ;

        class_<Bond>(m, "Bond")
            .def(init<SystemPtr, Id>())
            .def(self == self)
            .def(self != self)
            .def("__hash__", [](Bond const& a) { return hash(make_tuple(a.mol.get(), a.id)); })
            .def("__lt__", [](Bond& lhs, Bond& rhs) { return (lhs.mol < rhs.mol) || (lhs.id < rhs.id); })
            .def_readonly("_ptr", &Bond::mol)
            .def_readonly("id", &Bond::id, "unique id")
            .def("remove", [](Bond& b) { b.mol->delBond(b.id); }, "remove this Bond from the System")
            .def("otherId", [](Bond& b, Id i) { return b.mol->bond(b.id).other(i); })
            .def_property_readonly("i", [](Bond& b) { return b.mol->bond(b.id).i; }, "unique id of first atom")
            .def_property_readonly("j", [](Bond& b) { return b.mol->bond(b.id).j; }, "unique id of second atom")
            .def_property("order", [](Bond& b) { return b.mol->bond(b.id).order; },
                                   [](Bond& b, int o) { b.mol->bond(b.id).order=o; },
                                   "bond order (integer)")
            .def("__contains__", [](Bond& b, std::string const& key) { return BadId != b.mol->bondPropIndex(key); },
                    arg("key"), "does given custom Bond property exist?")
            .def("__getitem__", [](Bond& b, std::string const& key) {
                    Id col = b.mol->bondPropIndex(key);
                    if (bad(col)) {
                        PyErr_Format(PyExc_KeyError, "No such bond property '%s", key.data());
                        throw error_already_set();
                    }
                    return from_value_ref(b.mol->bondPropValue(b.id,col)); },
                    arg("key"), "get custom Bond property")
            .def("__setitem__", [](Bond& b, std::string const& key, object val) {
                    Id col = b.mol->bondPropIndex(key);
                    if (bad(col)) {
                        PyErr_Format(PyExc_KeyError, "No such bond property '%s", key.data());
                        throw error_already_set();
                    }
                    to_value_ref(val, b.mol->bondPropValue(b.id,col)); },
                    arg("key"), arg("val"), "set custom Bond property")
            ;

        class_<Residue>(m, "Residue")
            .def(init<SystemPtr, Id>())
            .def(self == self)
            .def(self != self)
            .def("__hash__", [](Residue const& a) { return hash(make_tuple(a.mol.get(), a.id)); })
            .def("__lt__", [](Residue& lhs, Residue& rhs) { return (lhs.mol < rhs.mol) || (lhs.id < rhs.id); })
            .def_readonly("_ptr", &Residue::mol)
            .def_readonly("id", &Residue::id, "unique id")
            .def_property("name", &Residue::name, &Residue::setName, "residue name")
            .def_property("insertion", &Residue::insertion, &Residue::setInsertion, "insertion code")
            .def_property("resid", &Residue::resid, &Residue::setResid, "PDB residue identifier")
            .def("remove", &Residue::remove, "remove this Residue from the System")
            .def("addAtom", &Residue::addAtom)
            .def("atoms", &Residue::atoms)
            .def_property_readonly("natoms", &Residue::natoms, "number of atoms in this residue")
            .def("chainId", &Residue::chainId, "id of parent Chain")
            .def("setChainId", &Residue::setChainId, "move residue to new Chain")
            ;

        class_<Chain>(m, "Chain")
            .def(init<SystemPtr, Id>())
            .def(self == self)
            .def(self != self)
            .def("__hash__", [](Chain const& a) { return hash(make_tuple(a.mol.get(), a.id)); })
            .def("__lt__", [](Chain& lhs, Chain& rhs) { return (lhs.mol < rhs.mol) || (lhs.id < rhs.id); })
            .def_readonly("_ptr", &Chain::mol)
            .def_readonly("id", &Chain::id, "unique id")
            .def_property("name", &Chain::name, &Chain::setName, "chain name")
            .def_property("segid", &Chain::segid, &Chain::setSegid, "segment name")
            .def_property_readonly("nresidues", &Chain::nresidues, "number of residues in chain")
            .def("ctId", &Chain::ctId, "id of parent Ct")
            .def("setCtId", &Chain::setCtId, "move Chain to new Ct")
            .def("remove", &Chain::remove, "remove this Chain from the Ct")
            .def("addResidue", &Chain::addResidue)
            .def("residues", &Chain::residues)
            ;

        class_<component_t>(m, "component_t")
            .def_property("name", &component_t::name, &component_t::setName)
            .def("keys", &component_t::keys)
            .def("get",  [](component_t& ct, String const& key) {
                    if (!ct.has(key)) { PyErr_SetString(PyExc_KeyError, key.data()); throw error_already_set(); }
                    return from_value_ref(ct.value(key));
                })
            .def("add",  [](component_t& ct, String const& key, object t) { ct.add(key, as_value_type(t)); })
            .def("set",  [](component_t& ct, String const& key, object val) { to_value_ref(val, ct.value(key)); })
            .def("type", [](component_t& ct, String const& key) -> handle { return ct.has(key) ? from_value_type(ct.type(key)) : none(); })
            .def("remove", &component_t::del)
            ;
    }
}}


