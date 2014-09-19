#include "unique_symbol.hxx"
#include "wrap_obj.hxx"
#include "schema.hxx"
#include "atomsel.hxx"
#include "append.hxx"
#include "clone.hxx"
#include "geom.hxx"
#include "contacts.hxx"

#include <numpy/ndarrayobject.h>

using namespace desres::msys;

namespace {

    atom_t& system_atom(System& m, Id id) { return m.atom(id); }
    bond_t& system_bond(System& m, Id id) { return m.bond(id); }
    residue_t& system_residue(System& m, Id id) { return m.residue(id); }
    chain_t& system_chain(System& m, Id id) { return m.chain(id); }
    component_t& system_ct(System& m, Id id) { return m.ct(id); }

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

    void destructor(PyObject* obj) { 
        Py_XDECREF(obj); 
    }
    typedef boost::shared_ptr<PyObject> objptr;

    PyObject* global_cell(object sysobj) {
        System& sys = extract<System&>(sysobj);
        npy_intp dims[2] = {3,3};
        PyObject* arr = PyArray_SimpleNewFromData(2,dims,NPY_FLOAT64,
                sys.global_cell[0]);
        Py_INCREF(PyArray_BASE(arr)=reinterpret_cast<PyObject*>(sysobj.ptr()));
        return arr;
    }

    PyObject* getpos3(object x, double **p) {
        PyObject* arr = PyArray_FromAny(
                x.ptr(), 
                PyArray_DescrFromType(NPY_FLOAT64),
                1, 1, NPY_C_CONTIGUOUS | NPY_ALIGNED,
                NULL);
        if (!arr) throw_error_already_set();
        if (PyArray_DIM(arr, 0) != 3) {
            PyErr_Format(PyExc_ValueError, "expected vector of length 3, got %ld", PyArray_DIM(arr, 0));
            throw_error_already_set();
        }
        if (p) *p = (double *)PyArray_DATA(arr);
        return arr;
    }

    double py_calc_distance(object A, object B) {
        double *a, *b;
        objptr pa(getpos3(A, &a), destructor);
        objptr pb(getpos3(B, &b), destructor);
        return calc_distance(a,b);
    }
    double py_calc_angle(object A, object B, object C) {
        double *a, *b, *c;
        objptr pa(getpos3(A, &a), destructor);
        objptr pb(getpos3(B, &b), destructor);
        objptr pc(getpos3(C, &c), destructor);
        return calc_angle(a,b,c);
    }
    double py_calc_dihedral(object A, object B, object C, object D) {
        double *a, *b, *c, *d;
        objptr pa(getpos3(A, &a), destructor);
        objptr pb(getpos3(B, &b), destructor);
        objptr pc(getpos3(C, &c), destructor);
        objptr pd(getpos3(D, &d), destructor);
        return calc_dihedral(a,b,c,d);
    }
    double py_calc_planarity(PyObject* obj) {
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
        return calc_planarity(n, pos);
    }

    bool py_line_intersects_tri(object A, object B, object C, 
                                object R, object S) {
        double *a, *b, *c, *r, *s;
        objptr pa(getpos3(A,&a), destructor);
        objptr pb(getpos3(B,&b), destructor);
        objptr pc(getpos3(C,&c), destructor);
        objptr pr(getpos3(R,&r), destructor);
        objptr ps(getpos3(S,&s), destructor);
        return line_intersects_tri(a,b,c,r,s);
    }

    void sys_translate(System& sys, double x, double y, double z) {
        for (Id i=0, n=sys.maxAtomId(); i<n; i++) {
            if (!sys.hasAtom(i)) continue;
            atom_t& atm = sys.atom(i);
            atm.x += x;
            atm.y += y;
            atm.z += z;
        }
    }

    list sys_topology(System& sys) {
        list result;
        for (Id i=0, n=sys.maxAtomId(); i<n; i++) {
            if (!sys.hasAtom(i)) continue;
            IdList atoms = sys.bondedAtoms(i);
            list sub;
            for (Id j=0, m=atoms.size(); j<m; j++) {
                sub.append(object(atoms[j]));
            }
            result.append(sub);
        }
        return result;
    }

    struct Finder {
        System const& mol;
        const bool skip_symmetric;
        list& contacts;
        Finder(System const& m, bool skip, list& c) 
        : mol(m), skip_symmetric(skip), contacts(c) {}
        bool exclude(Id i, Id j) const { 
            if (skip_symmetric && i>=j) return true;
            return !bad(mol.findBond(i,j)); 
        }
        void operator()(Id i, Id j, double d2) const {
            contacts.append(make_tuple(i,j,sqrt(d2)));
        }
    };

    list sys_find_contact_ids(System const& sys, double cutoff, 
                    object idobj, object otherobj, object posobj) {

        if (sys.atomCount()==0) return list();

        Id* idptr=NULL, *otherptr=NULL;
        unsigned nids=0, nother=0;
        const double* posptr=NULL;
        std::vector<double> posdata;
        IdList ids;
        objptr idarr, otherarr, posarr;
        bool skip_symmetric = true;

        if (idobj.ptr()==Py_None) {
            ids=sys.atoms();
            idptr = &ids[0];
            nids=ids.size();
        } else {
            idarr.reset(PyArray_FromAny(
                        idobj.ptr(),
                        PyArray_DescrFromType(NPY_UINT32),
                        1,1, NPY_C_CONTIGUOUS | NPY_FORCECAST, NULL),
                    destructor);
            if (!idarr) throw_error_already_set();
            idptr=(Id *)PyArray_DATA(idarr.get());
            nids=PyArray_DIM(idarr.get(), 0);
        }

        if (otherobj.ptr()==Py_None) {
            otherptr = idptr;
            nother=nids;
        } else {
            skip_symmetric = false;
            otherarr.reset(PyArray_FromAny(
                        otherobj.ptr(),
                        PyArray_DescrFromType(NPY_UINT32),
                        1,1, NPY_C_CONTIGUOUS | NPY_FORCECAST, NULL),
                    destructor);
            if (!otherarr) throw_error_already_set();
            otherptr = (Id *)PyArray_DATA(otherarr.get());
            nother = PyArray_DIM(otherarr.get(), 0);
        }

        if (posobj.ptr()==Py_None) {
            posdata.reserve(3*sys.atomCount());
            for (Id i=0, n=sys.maxAtomId(); i<n; i++) {
                if (!sys.hasAtom(i)) continue;
                atom_t const& atm = sys.atom(i);
                posdata.push_back(atm.x);
                posdata.push_back(atm.y);
                posdata.push_back(atm.z);
            }
            posptr= &posdata[0];
        } else {
            posarr.reset(PyArray_FromAny(
                        posobj.ptr(),
                        PyArray_DescrFromType(NPY_FLOAT64),
                        2,2,NPY_C_CONTIGUOUS, NULL),
                    destructor);
            if (!posarr) throw_error_already_set();
            posptr = (double *)PyArray_DATA(posarr.get());
        }

        list result;
        find_contacts(cutoff, posptr,
                      idptr, idptr+nids,
                      otherptr, otherptr+nother,
                      Finder(sys, skip_symmetric, result));
        return result;
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
        } else if (PyArray_Check(idobj.ptr())) {
            /* convert to int64 */
            PyObject* idarr = PyArray_FromAny(
                idobj.ptr(),
                PyArray_DescrFromType(NPY_INT64),
                1, 1,   /* must be 1-dimensional */
                NPY_C_CONTIGUOUS | NPY_ALIGNED,
                NULL);
            if (!idarr) throw_error_already_set();
            const int64_t* ids = (const int64_t *)PyArray_DATA(idarr);
            boost::shared_ptr<PyObject> _(idarr, destructor);

            /* allocate return array */
            dims[0] = PyArray_DIM(idarr,0);
            arr = PyArray_SimpleNew(2,dims,NPY_FLOAT64);
            if (!arr) throw_error_already_set();
            double* ptr = (double *)PyArray_DATA(arr);
            Py_ssize_t i,n = dims[0];
            for (i=0; i<n; i++) {
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

    PyObject* array_Atomselect(SystemPtr mol, std::string const& sel) {
        IdList ids = Atomselect(mol,sel);
        npy_intp dims[1];
        dims[0] = ids.size();
        PyObject *arr = PyArray_SimpleNew(1,dims,NPY_UINT32);
        if (!arr) throw_error_already_set();
        if (!ids.empty()) memcpy(PyArray_DATA(arr), &ids[0],
                ids.size()*sizeof(Id));
        return arr;
    }

    PyObject* list_atoms(SystemPtr mol) {
        Id i,n = mol->atomCount(), m = mol->maxAtomId();
        PyObject *L = PyList_New(n);
        if (!L) throw_error_already_set();
        if (n==m) for (i=0; i<m; i++) {
            PyList_SET_ITEM(L,i,PyInt_FromLong(i));
        } else    for (i=0, n=0; i<m; i++) {
            if (!mol->hasAtom(i)) continue;
            PyList_SET_ITEM(L,n++, PyInt_FromLong(i));
        }
        return L;
    }

    bool valid_permutation(SystemPtr mol, IdList ids) {
        /* check for duplicates */
        if (sort_unique(ids)) return false;
        /* check for dead atoms */
        BOOST_FOREACH(Id id, ids) if (!mol->hasAtom(id)) return false;
        /* check size */
        return ids.size() == mol->atomCount();
    }
}

namespace desres { namespace msys { 

    void export_system() {
        import_array();
        if (PyErr_Occurred()) return;
        
        def("calc_distance", py_calc_distance);
        def("calc_angle", py_calc_angle);
        def("calc_dihedral", py_calc_dihedral);
        def("calc_planarity", py_calc_planarity);
        def("line_intersects_tri", py_line_intersects_tri);

        class_<NonbondedInfo>("NonbondedInfo", no_init)
            .def_readwrite("vdw_funct", &NonbondedInfo::vdw_funct,
                    "Name of the vdw functional form; e.g., 'vdw_12_6'")
            .def_readwrite("vdw_rule", &NonbondedInfo::vdw_rule,
                    "Nonbonded combining rule; e.g., 'arithmetic/geometric'")
            .def_readwrite("es_funct", &NonbondedInfo::es_funct,
                    "Name of the electrostatic functional form")
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

        enum_<CloneOption::Flags>("CloneOption")
            .value("Default",       CloneOption::Default)
            .value("ShareParams",   CloneOption::ShareParams)
            ;

        def("Clone", Clone, 
                (arg("system"),
                 arg("ids"),
                 arg("flags")=CloneOption::Default));

        def("TableSchemas", table_schemas);
        def("NonbondedSchemas", nonbonded_schemas);

        class_<System,SystemPtr>("SystemPtr", no_init)
            .def("__eq__",      list_eq<SystemPtr>)
            .def("__ne__",      list_ne<SystemPtr>)
            .def("__hash__",    obj_hash<SystemPtr>)
            .def("create",  &System::create).staticmethod("create")
            .def_readwrite("name", &System::name)
            .add_property("global_cell", global_cell)
            .def_readwrite("nonbonded_info", &System::nonbonded_info)

            /* accessor */
            .def("atom",        system_atom, return_ptr())
            .def("bond",        system_bond, return_ptr())
            .def("residue",     system_residue, return_ptr())
            .def("chain",       system_chain, return_ptr())
            .def("ct",          system_ct, return_ptr())

            /* add element */
            .def("addAtom",     &System::addAtom)
            .def("addBond",     &System::addBond)
            .def("addResidue",  &System::addResidue)
            .def("addChain",    &System::addChain)
            .def("addCt",       &System::addCt)

            /* has element */
            .def("hasAtom",     &System::hasAtom)
            .def("hasBond",     &System::hasBond)
            .def("hasResidue",  &System::hasResidue)
            .def("hasChain",    &System::hasChain)
            .def("hasCt",       &System::hasCt)

            /* del element */
            .def("delAtom",     &System::delAtom)
            .def("delBond",     &System::delBond)
            .def("delResidue",  &System::delResidue)
            .def("delChain",    &System::delChain)
            .def("delCt",       &System::delCt)

            /* delete list of elements */
            .def("delAtoms",    &System::delAtoms<IdList>)
            .def("delBonds",    &System::delBonds<IdList>)
            .def("delResidues", &System::delResidues<IdList>)
            .def("delChains",   &System::delChains<IdList>)

            /* max element */
            .def("maxAtomId",   &System::maxAtomId)
            .def("maxBondId",   &System::maxBondId)
            .def("maxResidueId",&System::maxResidueId)
            .def("maxChainId",  &System::maxChainId)
            .def("maxCtId",     &System::maxCtId)

            /* list of element ids */
            .def("atoms",       &System::atoms)
            .def("atomsAsList", list_atoms)
            .def("bonds",       &System::bonds)
            .def("residues",    &System::residues)
            .def("chains",      &System::chains)
            .def("cts",         &System::cts)


            /* count of elements */
            .def("atomCount",   &System::atomCount)
            .def("bondCount",   &System::bondCount)
            .def("residueCount",&System::residueCount)
            .def("chainCount",  &System::chainCount)
            .def("ctCount",     &System::ctCount)

            /* count of subelements */
            .def("atomCountForResidue", &System::atomCountForResidue)
            .def("bondCountForAtom",    &System::bondCountForAtom)
            .def("residueCountForChain",&System::residueCountForChain)
            .def("chainCountForCt",     &System::chainCountForCt)
            .def("atomCountForCt",      &System::atomCountForCt)
            
            /* list of subelements */
            .def("atomsForResidue", &System::atomsForResidue, return_const() )
            .def("bondsForAtom",    &System::bondsForAtom, return_const())
            .def("residuesForChain",&System::residuesForChain, return_const())
            .def("chainsForCt",     &System::chainsForCt, return_const())
            .def("atomsForCt",      &System::atomsForCt)
            .def("bondsForCt",      &System::bondsForCt)

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

            /* chain props */
            .def("setCt",       &System::setCt)

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

            /* schemas */
            .def("addTableFromSchema", AddTable)
            .def("addNonbondedFromSchema", AddNonbonded)

            /* atom selection */
            .def("select", Atomselect)
            .def("selectAsList", list_Atomselect)
            .def("selectAsArray", array_Atomselect)

            /* append */
            .def("append", AppendSystem)

            /* miscellaneous */
            .def("validPermutation", valid_permutation)
            .def("orderedIds",    &System::orderedIds)
            .def("updateFragids", update_fragids)
            .def("analyze",     &System::analyze)
            .def("findBond",    &System::findBond)
            .def("provenance",      sys_provenance)
            .def("coalesceTables",    &System::coalesceTables)
            .def("translate",       sys_translate)
            .def("findContactIds",  sys_find_contact_ids)
            .def("topology", sys_topology)

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
