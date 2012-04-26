#include "wrap_obj.hxx"
#include "schema.hxx"
#include "atomsel.hxx"
#include "append.hxx"
#include "clone.hxx"
#include "mae.hxx"
#include "dms.hxx"
#include "pdb.hxx"
#include "amber.hxx"

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

    Provenance prov_from_args( object o ) {
        int argc=len(o);
        if (!argc) return Provenance();
        std::vector<std::string> strings;
        std::vector<char *> argv;
        for (int i=0; i<argc; i++) {
            strings.push_back(extract<std::string>(o[i]));
            argv.push_back(const_cast<char *>(strings.back().c_str()));
        }
        return Provenance::fromArgs(argc, &argv[0]);
    }

    list sys_provenance( System const& sys ) {
        list L;
        std::vector<Provenance> const& prov = sys.provenance();
        for (unsigned i=0; i<prov.size(); i++) {
            L.append(object(prov[i]));
        }
        return L;
    }

    SystemPtr import_mae_from_buffer(PyObject* obj, bool ignore_unrecognized) {
        Py_buffer view[1];
        if (PyObject_GetBuffer(obj, view, PyBUF_ND)) {
            throw_error_already_set();
        }
        boost::shared_ptr<Py_buffer> ptr(view, PyBuffer_Release);
        const char* bytes = reinterpret_cast<const char *>(view->buf);
        return ImportMAEFromBytes(bytes, view->len, ignore_unrecognized);
    }

    SystemPtr import_dms_from_buffer( PyObject* obj, bool structure_only ) {
        Py_buffer view[1];
        if (PyObject_GetBuffer(obj, view, PyBUF_ND)) {
            throw_error_already_set();
        }
        boost::shared_ptr<Py_buffer> ptr(view, PyBuffer_Release);
        char* bytes = reinterpret_cast<char *>(view->buf);
        return ImportDMSFromBytes(bytes, view->len, structure_only);
    }

    PyObject* sys_getpos(System const& sys) {
        PyObject *result = PyList_New(sys.atomCount());
        Id i,j=0,n = sys.maxAtomId();
        for (i=0; i<n; i++) {
            if (!sys.hasAtom(i)) continue;
            PyObject* xyz = PyList_New(3);
            PyList_SET_ITEM(xyz,0, PyFloat_FromDouble(sys.atom(i).x));
            PyList_SET_ITEM(xyz,1, PyFloat_FromDouble(sys.atom(i).y));
            PyList_SET_ITEM(xyz,2, PyFloat_FromDouble(sys.atom(i).z));
            PyList_SET_ITEM(result,j++,xyz);
        }
        return result;
    }

    void sys_setpos(System& sys, PyObject* obj) {
        PyObject* fast = PySequence_Fast(obj, "expected a sequence");
        if (!fast) throw_error_already_set();
        boost::shared_ptr<PyObject> ptr(fast, Py_DecRef);
        Id i, n = PySequence_Fast_GET_SIZE(fast);
        if (n!=sys.atomCount()) {
            PyErr_Format(PyExc_ValueError, "sequence has %d items, expected %d",
                    n, sys.atomCount());
            throw_error_already_set();
        }
        PyObject** items = PySequence_Fast_ITEMS(fast);
        for (i=0, n=0; i<sys.maxAtomId(); i++) {
            if (!sys.hasAtom(i)) continue;
            PyObject* item = PySequence_Fast(items[n++], 
                    "expected sequence of 3-tuples");
            if (!item) throw_error_already_set();
            boost::shared_ptr<PyObject> ptr(item, Py_DecRef);

            Id m = PySequence_Fast_GET_SIZE(item);
            if (m!=3) {
                PyErr_Format(PyExc_ValueError, 
                        "sequence item %d has %d items, expected 3", n, m);
                throw_error_already_set();
            }
            double x = PyFloat_AsDouble(PySequence_Fast_GET_ITEM(item,0));
            if (PyErr_Occurred()) throw_error_already_set();
            double y = PyFloat_AsDouble(PySequence_Fast_GET_ITEM(item,1));
            if (PyErr_Occurred()) throw_error_already_set();
            double z = PyFloat_AsDouble(PySequence_Fast_GET_ITEM(item,2));
            if (PyErr_Occurred()) throw_error_already_set();
            sys.atom(i).x=x;
            sys.atom(i).y=y;
            sys.atom(i).z=z;
        }
    }

    list glue_pairs(System const& sys) {
        list result;
        std::vector<glue_t> glue = sys.gluePairs();
        for (unsigned i=0; i<glue.size(); i++) {
            result.append(make_tuple(glue[i].first, glue[i].second));
        }
        return result;
    }

    /* FIXME: use a struct holding a member function implementing operator() */
    list table_names(System const& sys) {
        std::vector<std::string> names = sys.tableNames();
        list L;
        for (unsigned i=0; i<names.size(); i++) L.append(object(names[i]));
        return L;
    }
    list aux_table_names(System const& sys) {
        std::vector<std::string> names = sys.auxTableNames();
        list L;
        for (unsigned i=0; i<names.size(); i++) L.append(object(names[i]));
        return L;
    }
    list selection_macros(System const& sys) {
        std::vector<std::string> names = sys.selectionMacros();
        list L;
        for (unsigned i=0; i<names.size(); i++) L.append(object(names[i]));
        return L;
    }
    list table_schemas() {
        std::vector<std::string> names = TableSchemas();
        list L;
        for (unsigned i=0; i<names.size(); i++) L.append(object(names[i]));
        return L;
    }
    list nonbonded_schemas() {
        std::vector<std::string> names = NonbondedSchemas();
        list L;
        for (unsigned i=0; i<names.size(); i++) L.append(object(names[i]));
        return L;
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
                    "Name of the nonbonded functional form; e.g., \"vdw_12_6\"")
            .def_readwrite("vdw_rule", &NonbondedInfo::vdw_rule,
                    "Nonbonded combining rule; e.g., \"arithmetic/geometric\"")
            ;

        class_<Provenance>("Provenance", init<>())
            .def("fromArgs", prov_from_args).staticmethod("fromArgs")
            .def_readwrite("version", &Provenance::version)
            .def_readwrite("timestamp", &Provenance::timestamp)
            .def_readwrite("user", &Provenance::user)
            .def_readwrite("workdir", &Provenance::workdir)
            .def_readwrite("cmdline", &Provenance::cmdline)
            ;

        def("ImportDMS", ImportDMS);
        def("ImportDMSFromBuffer", import_dms_from_buffer);
        def("ExportDMS", ExportDMS);
        def("ImportMAE", ImportMAE);
        def("ImportMAEFromBuffer", import_mae_from_buffer);
        def("ExportMAE", ExportMAE);
        def("ImportPDB", ImportPDB);
        def("ImportPrmTop", ImportPrmTop);
        def("ImportCrdCoordinates", ImportCrdCoordinates);

        def("Clone", Clone);

        def("TableSchemas", table_schemas);
        def("NonbondedSchemas", nonbonded_schemas);

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
            .def("tableNames",  table_names)
            .def("tableName",   &System::tableName)
            .def("table",       &System::table)
            .def("addTable",    &System::addTable)
            .def("delTable",    &System::delTable)
            .def("removeTable", &System::removeTable)
            .def("renameTable", &System::renameTable)

            /* atom props */
            .def("setResidue",  &System::setResidue)
            .def("bondedAtoms", &System::bondedAtoms)

            /* residue props */
            .def("setChain",    &System::setChain)

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

            /* auxiliary tables */
            .def("auxTableNames",aux_table_names)
            .def("auxTable",     &System::auxTable)
            .def("addAuxTable",  &System::addAuxTable)
            .def("delAuxTable",  &System::delAuxTable)
            .def("removeAuxTable",&System::removeAuxTable)

            /* selection macros */
            .def("addSelectionMacro", &System::addSelectionMacro)
            .def("delSelectionMacro", &System::delSelectionMacro)
            .def("selectionMacroDefinition", &System::selectionMacroDefinition, return_const())
            .def("selectionMacroCount", &System::selectionMacroCount)
            .def("selectionMacros", selection_macros)

            /* schemas */
            .def("addTableFromSchema", AddTable)
            .def("addNonbondedFromSchema", AddNonbonded)

            /* atom selection */
            .def("select", Atomselect)

            /* append */
            .def("append", AppendSystem)

            /* glue */
            .def("glueCount", &System::glueCount)
            .def("gluePairs", glue_pairs)
            .def("hasGluePair", &System::hasGluePair)
            .def("addGluePair", &System::addGluePair)
            .def("delGluePair", &System::delGluePair)

            /* miscellaneous */
            .def("orderedIds",    &System::orderedIds)
            .def("updateFragids", update_fragids)
            .def("findBond",    &System::findBond)
            .def("provenance",      sys_provenance)
            .def("coalesceTables",    &System::coalesceTables)

            .def("getPositions",    sys_getpos)
            .def("setPositions",    sys_setpos)
            ;
    }

}}
