#include "wrap_obj.hxx"
#include "schema.hxx"
#include "atomsel.hxx"
#include "append.hxx"
#include "alchemical.hxx"
#include "clone.hxx"
#include "mae.hxx"
#include "dms.hxx"

#include <fstream>

using namespace desres::msys;

namespace {

    Vec3& get_A(GlobalCell& c) { return c.A; }
    Vec3& get_B(GlobalCell& c) { return c.B; }
    Vec3& get_C(GlobalCell& c) { return c.C; }
    Vec3& cell_getitem(GlobalCell& c, Py_ssize_t i) {
        if (i<0) i+=3;
        if (i<0 || i>=3) {
            PyErr_Format(PyExc_IndexError, "Index must be 0, 1 or 2");
            throw_error_already_set();
        }
        return c[i];
    }

    atom_t& system_atom(System& m, Id id) { return m.atom(id); }
    bond_t& system_bond(System& m, Id id) { return m.bond(id); }
    residue_t& system_residue(System& m, Id id) { return m.residue(id); }
    chain_t& system_chain(System& m, Id id) { return m.chain(id); }

    void del_atoms(System& sys, object elems) {
        list L(elems);
        IdList ids(len(L));
        for (unsigned i=0; i<ids.size(); i++) ids[i]=extract<Id>(L[i]);
        sys.delAtoms(ids);
    }
    void del_bonds(System& sys, object elems) {
        list L(elems);
        IdList ids(len(L));
        for (unsigned i=0; i<ids.size(); i++) ids[i]=extract<Id>(L[i]);
        sys.delBonds(ids);
    }
    void del_residues(System& sys, object elems) {
        list L(elems);
        IdList ids(len(L));
        for (unsigned i=0; i<ids.size(); i++) ids[i]=extract<Id>(L[i]);
        sys.delResidues(ids);
    }
    void del_chains(System& sys, object elems) {
        list L(elems);
        IdList ids(len(L));
        for (unsigned i=0; i<ids.size(); i++) ids[i]=extract<Id>(L[i]);
        sys.delChains(ids);
    }

    PyObject* atom_prop_type(System& sys, Id col) {
        return from_value_type(sys.atomPropType(col));
    }
    Id add_atom_prop( System& sys, String const& name, object type ) {
        return sys.addAtomProp(name, as_value_type(type));
    }
    object get_atom_prop(System& p, Id row, String const& name) {
        Id col = p.atomPropIndex(name);
        if (bad(col)) {
            PyErr_Format(PyExc_KeyError, 
                    "No such atom property '%s", name.c_str());
            throw_error_already_set();
        }
        return from_value_ref(p.atomPropValue(row,col));
    }
    void set_atom_prop(System& p, Id row, String const& name, object newval) {
        Id col = p.atomPropIndex(name);
        if (bad(col)) {
            PyErr_Format(PyExc_KeyError, 
                    "No such atom property '%s", name.c_str());
            throw_error_already_set();
        }
        to_value_ref(newval, p.atomPropValue(row,col));
    }

    PyObject* bond_prop_type(System& sys, Id col) {
        return from_value_type(sys.bondPropType(col));
    }
    Id add_bond_prop( System& sys, String const& name, object type ) {
        return sys.addBondProp(name, as_value_type(type));
    }
    object get_bond_prop(System& p, Id row, String const& name) {
        Id col = p.bondPropIndex(name);
        if (bad(col)) {
            PyErr_Format(PyExc_KeyError, 
                    "No such bond property '%s", name.c_str());
            throw_error_already_set();
        }
        return from_value_ref(p.bondPropValue(row,col));
    }
    void set_bond_prop(System& p, Id row, String const& name, object newval) {
        Id col = p.bondPropIndex(name);
        if (bad(col)) {
            PyErr_Format(PyExc_KeyError, 
                    "No such bond property '%s", name.c_str());
            throw_error_already_set();
        }
        to_value_ref(newval, p.bondPropValue(row,col));
    }


    list update_fragids(System& p) {
        MultiIdList fragments;
        p.updateFragids(&fragments);
        list result;
        for (unsigned i=0; i<fragments.size(); i++) {
            result.append(object(fragments[i]));
        }
        return result;
    }
}

namespace desres { namespace msys { 

    void export_system() {
        
        class_<GlobalCell>("GlobalCell", no_init)
            .add_property("A", make_function(get_A, return_ptr()))
            .add_property("B", make_function(get_B, return_ptr()))
            .add_property("C", make_function(get_C, return_ptr()))
            .def("__getitem__", cell_getitem, return_ptr())
            ;

        class_<NonbondedInfo>("NonbondedInfo", no_init)
            .def_readwrite("vdw_funct", &NonbondedInfo::vdw_funct,
                    "Name of the nonbonded functional form; e.g., vdw_12_6")
            .def_readwrite("vdw_rule", &NonbondedInfo::vdw_rule,
                    "Nonbonded combining rule; e.g., arithmetic/geometric")
            ;

        def("ImportDMS", ImportDMS);
        def("ExportDMS", ExportDMS);
        def("ImportMAE", ImportMAE);

        def("Clone", Clone);
        def("CreateAlchemical", CreateAlchemical);

        class_<System,SystemPtr>("SystemPtr", no_init)
            .def("__eq__",      list_eq<SystemPtr>)
            .def("__ne__",      list_ne<SystemPtr>)
            .def("__hash__",    obj_hash<SystemPtr>)
            .def("create",  &System::create).staticmethod("create")
            .def_readwrite("name", &System::name)
            .def_readwrite("global_cell", &System::global_cell)
            .def_readwrite("nonbonded_info", &System::nonbonded_info)

            /* accessor */
            .def("atom",        system_atom, return_ptr())
            .def("bond",        system_bond, return_ptr())
            .def("residue",     system_residue, return_ptr())
            .def("chain",       system_chain, return_ptr())

            /* add element */
            .def("addAtom",     &System::addAtom)
            .def("addBond",     &System::addBond)
            .def("addResidue",  &System::addResidue)
            .def("addChain",    &System::addChain)

            /* has element */
            .def("hasAtom",     &System::hasAtom)
            .def("hasBond",     &System::hasBond)
            .def("hasResidue",  &System::hasResidue)
            .def("hasChain",    &System::hasChain)

            /* del element */
            .def("delAtom",     &System::delAtom)
            .def("delBond",     &System::delBond)
            .def("delResidue",  &System::delResidue)
            .def("delChain",    &System::delChain)

            /* delete list of elements */
            .def("delAtoms",    del_atoms)
            .def("delBonds",    del_bonds)
            .def("delResidues", del_residues)
            .def("delChains",   del_chains)

            /* max element */
            .def("maxAtomId",   &System::maxAtomId)
            .def("maxBondId",   &System::maxBondId)
            .def("maxResidueId",&System::maxResidueId)
            .def("maxChainId",  &System::maxChainId)

            /* list of element ids */
            .def("atoms",       &System::atoms)
            .def("bonds",       &System::bonds)
            .def("residues",    &System::residues)
            .def("chains",      &System::chains)

            /* count of elements */
            .def("atomCount",   &System::atomCount)
            .def("bondCount",   &System::bondCount)
            .def("residueCount",&System::residueCount)
            .def("chainCount",  &System::chainCount)

            /* count of subelements */
            .def("atomCountForResidue", &System::atomCountForResidue)
            .def("bondCountForAtom",    &System::bondCountForAtom)
            .def("residueCountForChain",&System::residueCountForChain)
            
            /* list of subelements */
            .def("atomsForResidue", &System::atomsForResidue, return_const() )
            .def("bondsForAtom",    &System::bondsForAtom, return_const())
            .def("residuesForChain",&System::residuesForChain, return_const())

            /* tables */
            .def("tableNames",  &System::tableNames)
            .def("tableName",   &System::tableName)
            .def("table",       &System::table)
            .def("addTable",    &System::addTable)
            .def("delTable",    &System::delTable)
            .def("removeTable", &System::removeTable)
            .def("renameTable", &System::renameTable)

            /* atom props */
            .def("setResidue",  &System::setResidue)
            .def("bondedAtoms", &System::bondedAtoms)

            /* extended atom props */
            .def("atomPropCount",&System::atomPropCount)
            .def("atomPropName", &System::atomPropName)
            .def("atomPropIndex",&System::atomPropIndex)
            .def("atomPropType", atom_prop_type)
            .def("addAtomProp",  add_atom_prop)
            .def("delAtomProp",  &System::delAtomProp)
            .def("getAtomProp",  get_atom_prop)
            .def("setAtomProp",  set_atom_prop)

            /* extended bond props */
            .def("bondPropCount",&System::bondPropCount)
            .def("bondPropName", &System::bondPropName)
            .def("bondPropIndex",&System::bondPropIndex)
            .def("bondPropType", bond_prop_type)
            .def("addBondProp",  add_bond_prop)
            .def("delBondProp",  &System::delBondProp)
            .def("getBondProp",  get_bond_prop)
            .def("setBondProp",  set_bond_prop)

            /* extras */
            .def("extraNames",  &System::extraNames)
            .def("extra",       &System::extra)
            .def("addExtra",    &System::addExtra)
            .def("delExtra",    &System::delExtra)
            .def("removeExtra", &System::removeExtra)

            /* schemas */
            .def("addTableFromSchema", AddTable)
            .def("addNonbondedFromSchema", AddNonbonded)

            /* atom selection */
            .def("atomselect", Atomselect)

            /* append */
            .def("append", AppendSystem)

            /* miscellaneous */
            .def("reassignGids",    &System::reassignGids)
            .def("updateFragids", update_fragids)
            .def("findBond",    &System::findBond)
            ;
    }

}}
