#include "wrap_obj.hxx"
#include "schema.hxx"
#include "atomsel.hxx"
#include "append.hxx"
#include "clone.hxx"
#include "mae.hxx"
#include "dms.hxx"
#include "pdb.hxx"
#include "load.hxx"
#include "amber.hxx"

#include <fstream>
#include <numpy/ndarrayobject.h>

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

    PyObject* sys_getpos(System const& sys, object idobj) {
        npy_intp dims[2];
        dims[1] = 3;
        PyObject* arr=NULL;
        if (idobj.ptr()==Py_None) {
            dims[0] = sys.atomCount();
            arr = PyArray_SimpleNew(2,dims,NPY_FLOAT64);
            if (!arr) throw_error_already_set();
            double* ptr = (double *)PyArray_DATA(arr);
            for (Id i=0, n=sys.maxAtomId(); i<n; i++) {
                if (!sys.hasAtom(i)) continue;
                atom_t const& atm = sys.atom(i);
                ptr[0] = atm.x;
                ptr[1] = atm.y;
                ptr[2] = atm.z;
                ptr += 3;
            }
        } else {
            IdList const& ids = extract<IdList const&>(idobj);
            dims[0] = ids.size();
            arr = PyArray_SimpleNew(2,dims,NPY_FLOAT64);
            if (!arr) throw_error_already_set();
            double* ptr = (double *)PyArray_DATA(arr);
            for (Id i=0, n = ids.size(); i<n; i++) {
                Id id = ids[i];
                if (!sys.hasAtom(id)) {
                    Py_DECREF(arr);
                    PyErr_Format(PyExc_ValueError, "Invalid id %u", id);
                    throw_error_already_set();
                }
                atom_t const& atm = sys.atom(id);
                ptr[0] = atm.x;
                ptr[1] = atm.y;
                ptr[2] = atm.z;
                ptr += 3;
            }
        }
        return arr;
    }

    PyObject* sys_getvel(System const& sys, object idobj) {
        npy_intp dims[2];
        dims[1] = 3;
        PyObject* arr=NULL;
        if (idobj.ptr()==Py_None) {
            dims[0] = sys.atomCount();
            arr = PyArray_SimpleNew(2,dims,NPY_FLOAT64);
            if (!arr) throw_error_already_set();
            double* ptr = (double *)PyArray_DATA(arr);
            for (Id i=0, n=sys.maxAtomId(); i<n; i++) {
                if (!sys.hasAtom(i)) continue;
                atom_t const& atm = sys.atom(i);
                ptr[0] = atm.vx;
                ptr[1] = atm.vy;
                ptr[2] = atm.vz;
                ptr += 3;
            }
        } else {
            IdList const& ids = extract<IdList const&>(idobj);
            dims[0] = ids.size();
            arr = PyArray_SimpleNew(2,dims,NPY_FLOAT64);
            if (!arr) throw_error_already_set();
            double* ptr = (double *)PyArray_DATA(arr);
            for (Id i=0, n = ids.size(); i<n; i++) {
                Id id = ids[i];
                if (!sys.hasAtom(id)) {
                    Py_DECREF(arr);
                    PyErr_Format(PyExc_ValueError, "Invalid id %u", id);
                    throw_error_already_set();
                }
                atom_t const& atm = sys.atom(id);
                ptr[0] = atm.vx;
                ptr[1] = atm.vy;
                ptr[2] = atm.vz;
                ptr += 3;
            }
        }
        return arr;
    }

    void destructor(PyObject* obj) { 
        Py_DECREF(obj); 
    }

    void sys_setpos(System& sys, PyObject* obj, object idobj) {
        PyObject* arr = PyArray_FromAny(
                obj,
                PyArray_DescrFromType(NPY_FLOAT64),
                2, 2,   /* must be 2-dimensional */
                NPY_C_CONTIGUOUS | NPY_ALIGNED,
                NULL);
        if (!arr) throw_error_already_set();
        boost::shared_ptr<PyObject> _(arr, destructor);
        Id n = PyArray_DIM(arr,0);
        if (PyArray_DIM(arr,1)!=3) {
            PyErr_Format(PyExc_ValueError, 
                    "Supplied %ld-d positions, expected 3-d", 
                    PyArray_DIM(arr,1));
            throw_error_already_set();
        }
        const double* pos = (const double *)PyArray_DATA(arr);

        if (idobj.ptr()==Py_None) {
            if (n!=sys.atomCount()) {
                PyErr_Format(PyExc_ValueError, 
                        "Supplied %u positions, but system has %u atoms", 
                        n, sys.atomCount());
                throw_error_already_set();
            }
            for (Id i=0, n=sys.maxAtomId(); i<n; i++) {
                if (!sys.hasAtom(i)) continue;
                atom_t& atm = sys.atom(i);
                atm.x = pos[0];
                atm.y = pos[1];
                atm.z = pos[2];
                pos += 3;
            }
        } else {
            IdList const& ids = extract<IdList const&>(idobj);
            if (n != ids.size()) {
                PyErr_Format(PyExc_ValueError,
                        "Supplied %u positions != %lu ids", n, ids.size());
                throw_error_already_set();
            }
            for (Py_ssize_t i=0; i<n; i++) {
                Id id = ids[i];
                if (!sys.hasAtom(id)) {
                    PyErr_Format(PyExc_ValueError,
                            "Id %u at position %ld is invalid", id, i);
                    throw_error_already_set();
                }
                atom_t& atm = sys.atom(id);
                atm.x = pos[0];
                atm.y = pos[1];
                atm.z = pos[2];
                pos += 3;
            }
        }
    }

    void sys_setvel(System& sys, PyObject* obj, object idobj) {
        PyObject* arr = PyArray_FromAny(
                obj,
                PyArray_DescrFromType(NPY_FLOAT64),
                2, 2,   /* must be 2-dimensional */
                NPY_C_CONTIGUOUS | NPY_ALIGNED,
                NULL);
        if (!arr) throw_error_already_set();
        boost::shared_ptr<PyObject> _(arr, destructor);
        Id n = PyArray_DIM(arr,0);
        if (PyArray_DIM(arr,1)!=3) {
            PyErr_Format(PyExc_ValueError, 
                    "Supplied %ld-d velocities, expected 3-d", 
                    PyArray_DIM(arr,1));
            throw_error_already_set();
        }
        const double* pos = (const double *)PyArray_DATA(arr);

        if (idobj.ptr()==Py_None) {
            if (n!=sys.atomCount()) {
                PyErr_Format(PyExc_ValueError, 
                        "Supplied %u velocities, but system has %u atoms", 
                        n, sys.atomCount());
                throw_error_already_set();
            }
            for (Id i=0, n=sys.maxAtomId(); i<n; i++) {
                if (!sys.hasAtom(i)) continue;
                atom_t& atm = sys.atom(i);
                atm.vx = pos[0];
                atm.vy = pos[1];
                atm.vz = pos[2];
                pos += 3;
            }
        } else {
            IdList const& ids = extract<IdList const&>(idobj);
            if (n != ids.size()) {
                PyErr_Format(PyExc_ValueError,
                        "Supplied %u velocities != %lu ids", n, ids.size());
                throw_error_already_set();
            }
            for (Py_ssize_t i=0; i<n; i++) {
                Id id = ids[i];
                if (!sys.hasAtom(id)) {
                    PyErr_Format(PyExc_ValueError,
                            "Id %u at position %ld is invalid", id, i);
                    throw_error_already_set();
                }
                atom_t& atm = sys.atom(id);
                atm.vx = pos[0];
                atm.vy = pos[1];
                atm.vz = pos[2];
                pos += 3;
            }
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

    PyObject* list_Atomselect(SystemPtr mol, std::string const& sel) {
        IdList ids = Atomselect(mol,sel);
        PyObject *L = PyList_New(ids.size());
        if (!L) throw_error_already_set();
        for (unsigned i=0; i<ids.size(); i++) {
            PyList_SET_ITEM(L,i,PyInt_FromLong(ids[i]));
        }
        return L;
    }
}

namespace desres { namespace msys { 

    void export_system() {
        import_array();
        if (PyErr_Occurred()) return;
        
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
            .def_readwrite("executable", &Provenance::executable)
            ;

        def("ImportDMS", ImportDMS);
        def("ImportDMSFromBuffer", import_dms_from_buffer);
        def("ExportDMS", ExportDMS);
        def("ImportMAE", ImportMAE);
        def("ImportMAEFromBuffer", import_mae_from_buffer);
        def("ExportMAE", ExportMAE);
        def("ImportPDB", ImportPDB);
        def("ExportPDB", ExportPDB);
        def("ImportPrmTop", ImportPrmTop);
        def("ImportCrdCoordinates", ImportCrdCoordinates);
        def("Load", Load,
                (arg("path"),
                 arg("opt_format")=object()));

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
            .def("selectAsList", list_Atomselect)

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
            .def("analyze",     &System::analyze)
            .def("findBond",    &System::findBond)
            .def("provenance",      sys_provenance)
            .def("coalesceTables",    &System::coalesceTables)

            .def("getPositions", sys_getpos,
                    (arg("ids")=object()))
            .def("setPositions",    sys_setpos,
                    (arg("pos"),
                     arg("ids")=object()))
            .def("getVelocities", sys_getvel,
                    (arg("ids")=object()))
            .def("setVelocities",    sys_setvel,
                    (arg("vel"),
                     arg("ids")=object()))
            ;
    }

}}
