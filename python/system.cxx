#include "unique_symbol.hxx"
#include "wrap_obj.hxx"
#include "capsule.hxx"

#include "schema.hxx"
#include "atomsel.hxx"
#include "append.hxx"
#include "clone.hxx"
#include "geom.hxx"
#include "dms.hxx"
#include "hbond.hxx"
#include "contacts.hxx"

#include <msys/hash.hxx>

#include <numpy/ndarrayobject.h>
#include <pfx/graph.hxx>
#include <pfx/cell.hxx>

using namespace desres::msys;

namespace {
#if PY_MAJOR_VERSION >= 3
    auto py_from_long = PyLong_FromLong;
    auto py_as_bytes = PyBytes_FromStringAndSize;
#else
    auto py_from_long = PyInt_FromLong;
    auto py_as_bytes = PyString_FromStringAndSize;
#endif
}

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
        return to_python(fragments);
    }

    Provenance prov_from_args( object o ) {
        int argc=len(o);
        if (!argc) return Provenance();
        std::vector<std::string> strings;
        std::vector<char *> argv;
        for (int i=0; i<argc; i++) {
            strings.push_back(extract<std::string>(o[i]));
        }
        for (int i=0; i<argc; i++) {
            argv.push_back(const_cast<char *>(strings[i].c_str()));
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
    typedef std::shared_ptr<PyObject> objptr;

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

    PyObject* make_projection(object cell) {
        npy_intp dims[2] = {3,3};
        PyObject* proj = PyArray_SimpleNew(2, dims, NPY_FLOAT64);
        pfx::make_projection((double*)PyArray_DATA(cell.ptr()),
                             (double*)PyArray_DATA(proj));
        return proj;
    }

    object wrap_vector(object cell, object proj, object pos) {
        pfx::wrap_vector<double>((double*)PyArray_DATA(cell.ptr()),
                                 (double*)PyArray_DATA(proj.ptr()),
                                 (double*)PyArray_DATA(pos.ptr()));
        return pos;
    }

    object wrap_vector_array(object cell, object proj, object pos) {
        pfx::wrap_vector_array<double>((double*)PyArray_DATA(cell.ptr()),
                                       (double*)PyArray_DATA(proj.ptr()),
                                       (double*)PyArray_DATA(pos.ptr()),
                                       extract<unsigned>(pos.attr("size")) / 3
                                       );
        return pos;
    }


    double py_calc_distance(object A, object B) {
        double *a, *b;
        objptr pa(getpos3(A, &a), destructor);
        objptr pb(getpos3(B, &b), destructor);
        return calc_distance(a,b);
    }
    double py_calc_vec_angle(object A, object B) {
        double *a, *b;
        objptr pa(getpos3(A, &a), destructor);
        objptr pb(getpos3(B, &b), destructor);
        return calc_vec_angle(a,b);
    }
    double py_calc_angle(object A, object B, object C) {
        double *a, *b, *c;
        objptr pa(getpos3(A, &a), destructor);
        objptr pb(getpos3(B, &b), destructor);
        objptr pc(getpos3(C, &c), destructor);
        return calc_angle(a,b,c);
    }
    double py_calc_vec_dihedral(object A, object B, object C) {
        double *a, *b, *c;
        objptr pa(getpos3(A, &a), destructor);
        objptr pb(getpos3(B, &b), destructor);
        objptr pc(getpos3(C, &c), destructor);
        return calc_vec_dihedral(a,b,c);
    }
    double py_calc_dihedral(object A, object B, object C, object D) {
        double *a, *b, *c, *d;
        objptr pa(getpos3(A, &a), destructor);
        objptr pb(getpos3(B, &b), destructor);
        objptr pc(getpos3(C, &c), destructor);
        objptr pd(getpos3(D, &d), destructor);
        return calc_dihedral(a,b,c,d);
    }

    PyObject* py_apply_dihedral_geometry(object A, object B, object C,
                                      double r, double theta, double phi) {
        double *a, *b, *c;
        objptr pa(getpos3(A, &a), destructor);
        objptr pb(getpos3(B, &b), destructor);
        objptr pc(getpos3(C, &c), destructor);
        npy_intp dim = 3;
        PyObject* arr = PyArray_SimpleNew(1,&dim,NPY_FLOAT64);
        double *d = (double *)PyArray_DATA(arr);
        apply_dihedral_geometry(d,a,b,c,r,theta,phi);
        return arr;
    }

    double py_calc_planarity(PyObject* obj) {
        PyObject* arr = PyArray_FromAny(
                obj,
                PyArray_DescrFromType(NPY_FLOAT64),
                2, 2,   /* must be 2-dimensional */
                NPY_C_CONTIGUOUS | NPY_ALIGNED,
                NULL);
        if (!arr) throw_error_already_set();
        std::shared_ptr<PyObject> _(arr, destructor);
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
        } else {
            /* convert to int64 */
            PyObject* idarr = PyArray_FromAny(
                idobj.ptr(),
                PyArray_DescrFromType(NPY_INT64),
                1, 1,   /* must be 1-dimensional */
                NPY_C_CONTIGUOUS | NPY_ALIGNED,
                NULL);
            if (!idarr) throw_error_already_set();
            const int64_t* ids = (const int64_t *)PyArray_DATA(idarr);
            std::shared_ptr<PyObject> _(idarr, destructor);

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
            /* convert to int64 */
            PyObject* idarr = PyArray_FromAny(
                idobj.ptr(),
                PyArray_DescrFromType(NPY_INT64),
                1, 1,   /* must be 1-dimensional */
                NPY_C_CONTIGUOUS | NPY_ALIGNED,
                NULL);
            if (!idarr) throw_error_already_set();
            const int64_t* ids = (const int64_t *)PyArray_DATA(idarr);
            std::shared_ptr<PyObject> _(idarr, destructor);

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
                ptr[0] = atm.vx;
                ptr[1] = atm.vy;
                ptr[2] = atm.vz;
                ptr += 3;
            }
        }
        return arr;
    }

    void sys_setpos(System& sys, PyObject* obj, PyObject* idobj) {
        PyObject* arr = PyArray_FromAny(
                obj,
                PyArray_DescrFromType(NPY_FLOAT64),
                2, 2,   /* must be 2-dimensional */
                NPY_C_CONTIGUOUS | NPY_ALIGNED,
                NULL);
        if (!arr) throw_error_already_set();
        std::shared_ptr<PyObject> _(arr, destructor);
        Id n = PyArray_DIM(arr,0);
        if (PyArray_DIM(arr,1)!=3) {
            PyErr_Format(PyExc_ValueError, 
                    "Supplied %ld-d positions, expected 3-d", 
                    PyArray_DIM(arr,1));
            throw_error_already_set();
        }
        const double* pos = (const double *)PyArray_DATA(arr);

        if (idobj==Py_None) {
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
            PyObject* idarr = PyArray_FromAny(
                    idobj,
                    PyArray_DescrFromType(NPY_INT64),
                    1,1,
                    NPY_C_CONTIGUOUS | NPY_ALIGNED,
                    NULL);
            if (!idarr) throw_error_already_set();
            std::shared_ptr<PyObject> _(idarr, destructor);
            if (n != (Id)PyArray_DIM(idarr,0)) {
                PyErr_Format(PyExc_ValueError,
                        "Supplied %u positions != %u ids",
                        n, (Id)PyArray_DIM(idarr,0));
                throw_error_already_set();
            }
            const Id* ids = (const Id*)PyArray_DATA(idarr);
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

    void sys_setvel(System& sys, PyObject* obj, PyObject* idobj) {
        PyObject* arr = PyArray_FromAny(
                obj,
                PyArray_DescrFromType(NPY_FLOAT64),
                2, 2,   /* must be 2-dimensional */
                NPY_C_CONTIGUOUS | NPY_ALIGNED,
                NULL);
        if (!arr) throw_error_already_set();
        std::shared_ptr<PyObject> _(arr, destructor);
        Id n = PyArray_DIM(arr,0);
        if (PyArray_DIM(arr,1)!=3) {
            PyErr_Format(PyExc_ValueError, 
                    "Supplied %ld-d velocities, expected 3-d", 
                    PyArray_DIM(arr,1));
            throw_error_already_set();
        }
        const double* pos = (const double *)PyArray_DATA(arr);

        if (idobj==Py_None) {
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
            PyObject* idarr = PyArray_FromAny(
                    idobj,
                    PyArray_DescrFromType(NPY_INT64),
                    1,1,
                    NPY_C_CONTIGUOUS | NPY_ALIGNED,
                    NULL);
            if (!idarr) throw_error_already_set();
            std::shared_ptr<PyObject> _(idarr, destructor);
            if (n != (Id)PyArray_DIM(idarr,0)) {
                PyErr_Format(PyExc_ValueError,
                        "Supplied %u velocities != %u ids",
                        n, (Id)PyArray_DIM(idarr,0));
                throw_error_already_set();
            }
            const Id* ids = (const Id*)PyArray_DATA(idarr);
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

    IdList wrap_atomselect(SystemPtr mol, std::string const& sel,
                           PyObject* posobj, PyObject* boxobj) {
        objptr posarr, boxarr;
        float* pos = NULL;
        double* box = NULL;
        if (posobj != Py_None) {
            posarr.reset(PyArray_FromAny(
                        posobj, PyArray_DescrFromType(NPY_FLOAT32),
                        2,2, NPY_C_CONTIGUOUS, NULL), destructor);
            if (!posarr) throw_error_already_set();
            if (PyArray_DIM(posarr.get(),0)!=mol->atomCount() ||
                PyArray_DIM(posarr.get(),1)!=3) {
                PyErr_Format(PyExc_ValueError, "pos has wrong shape");
                throw_error_already_set();
            }
            pos = static_cast<float*>(PyArray_DATA(posarr.get()));
        }
        if (boxobj != Py_None) {
            boxarr.reset(PyArray_FromAny(
                        boxobj, PyArray_DescrFromType(NPY_FLOAT64),
                        2,2, NPY_C_CONTIGUOUS, NULL), destructor);
            if (!boxarr) throw_error_already_set();
            if (PyArray_DIM(boxarr.get(),0)!=3 ||
                PyArray_DIM(boxarr.get(),1)!=3) {
                PyErr_Format(PyExc_ValueError, "box has wrong shape");
                throw_error_already_set();
            }
            box = static_cast<double*>(PyArray_DATA(boxarr.get()));
        }
        return Atomselect(mol, sel, pos, box);
    }

    PyObject* list_Atomselect(SystemPtr mol, std::string const& sel,
                               PyObject* pos, PyObject* box) {
        IdList ids = wrap_atomselect(mol,sel, pos, box);
        PyObject *L = PyList_New(ids.size());
        if (!L) throw_error_already_set();
        for (unsigned i=0; i<ids.size(); i++) {
            PyList_SET_ITEM(L,i,py_from_long(ids[i]));
        }
        return L;
    }

    PyObject* array_Atomselect(SystemPtr mol, std::string const& sel,
                               PyObject* pos, PyObject* box) {
        IdList ids = wrap_atomselect(mol,sel, pos, box);
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
            PyList_SET_ITEM(L,i,py_from_long(i));
        } else    for (i=0, n=0; i<m; i++) {
            if (!mol->hasAtom(i)) continue;
            PyList_SET_ITEM(L,n++, py_from_long(i));
        }
        return L;
    }

    list list_bonds(System& mol) { return to_python(mol.bonds()); }
    list list_residues(System& mol) { return to_python(mol.residues()); }
    list list_chains(System& mol) { return to_python(mol.chains()); }
    list list_cts(System& mol) { return to_python(mol.cts()); }

    list bonded_atoms(System& mol,Id i) { return to_python(mol.bondedAtoms(i));}
    list residue_atoms(System& mol,Id i) { return to_python(mol.atomsForResidue(i)); }
    list atom_bonds(System& mol,Id i) { return to_python(mol.bondsForAtom(i)); }
    list chain_residues(System& mol,Id i) { return to_python(mol.residuesForChain(i)); }
    list ct_chains(System& mol,Id i) { return to_python(mol.chainsForCt(i)); }
    list ct_atoms(System& mol,Id i) { return to_python(mol.atomsForCt(i)); }
    list ct_bonds(System& mol,Id i) { return to_python(mol.bondsForCt(i)); }

    objptr get_vec3d(PyObject* obj) {
        if (obj==Py_None) return objptr();
        PyObject* arr = PyArray_FromAny(
                obj,
                PyArray_DescrFromType(NPY_FLOAT64),
                1,1,
                NPY_C_CONTIGUOUS,
                NULL);
        if (!arr) throw_error_already_set();
        objptr ptr(arr, destructor);
        if (PyArray_DIM(arr,0)!=3) {
            PyErr_Format(PyExc_ValueError,
                    "Expected 3 elements in vector, got %ld",
                    PyArray_DIM(arr,0));
            throw_error_already_set();
        }
        return ptr;
    }

    HydrogenBond *init_hbond(
            PyObject* dobj,
            PyObject* aobj,
            PyObject* hobj,
            PyObject* cobj,
            PyObject* caobj) {

        objptr darr = get_vec3d(dobj);
        objptr aarr = get_vec3d(aobj);
        objptr harr = get_vec3d(hobj);
        objptr carr = get_vec3d(cobj);
        objptr caarr = get_vec3d(caobj);
        return new HydrogenBond(
                darr ? (const double *)PyArray_DATA(darr.get()) : NULL,
                aarr ? (const double *)PyArray_DATA(aarr.get()) : NULL,
                harr ? (const double *)PyArray_DATA(harr.get()) : NULL,
                carr ? (const double *)PyArray_DATA(carr.get()) : NULL,
                caarr? (const double *)PyArray_DATA(caarr.get()): NULL);
    }

    pfx::Graph* sys_topology(SystemPtr mol) { return new pfx::Graph(mol); }

    TermTablePtr wrap_system_add_table(SystemPtr mol, std::string const& name,
            Id natoms, PyObject* obj) {
        ParamTablePtr params;
        if (obj!=Py_None) {
            params = extract<ParamTablePtr>(obj);
        }
        return mol->addTable(name, natoms, params);
    }
    list ordered_ids(System& mol) {
        return to_python(mol.orderedIds());
    }

    SystemPtr wrap_clone(SystemPtr mol, list ids, CloneOption::Flags flags) {
        return Clone(mol, ids_from_python(ids), flags);
    }
    list append_system(SystemPtr dst, SystemPtr src, Id ct=BadId) {
        return to_python(AppendSystem(dst, src, ct));
    }

    struct system_pickle_suite : pickle_suite {
        static tuple getinitargs(SystemPtr mol) {
            std::string contents = FormatDMS(mol, Provenance());
            auto bytes = py_as_bytes(contents.data(), contents.size());
            return boost::python::make_tuple(
                    "dmscontents", handle<>(bytes));
        }
        static tuple getstate(SystemPtr mol) {
            return tuple();
        }
        static void setstate(SystemPtr mol, tuple s) {
        }
    };

    SystemPtr init_from_pickle(std::string const& format, std::string const& contents) {
        if (format=="dmscontents") {
            return ImportDMSFromBytes(contents.data(), contents.size());
        }
        throw std::runtime_error("Unsupported pickle format " + format);
    }
}

namespace desres { namespace msys { 

    void export_system() {
        _import_array();
        if (PyErr_Occurred()) return;
        
        def("bad", bad);
        scope().attr("BadId") = (Id)BadId;

        def("make_projection", make_projection);
        def("wrap_vector", wrap_vector);
        def("wrap_vector_array", wrap_vector_array);
        def("calc_distance", py_calc_distance);
        def("calc_vec_angle", py_calc_vec_angle);
        def("calc_angle", py_calc_angle);
        def("calc_vec_dihedral", py_calc_vec_dihedral);
        def("calc_dihedral", py_calc_dihedral);
        def("calc_planarity", py_calc_planarity);
        def("line_intersects_tri", py_line_intersects_tri);
        def("apply_dihedral_geometry", py_apply_dihedral_geometry);

        class_<pfx::Graph>("Topology", init<unsigned>())
            .add_property("nverts", &pfx::Graph::nverts)
            .add_property("nedges", &pfx::Graph::nedges)
            ;

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

            /* pickle support */
            .def_pickle(system_pickle_suite())
            .def("__init__", make_constructor(init_from_pickle))

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

            /* max element */
            .def("maxAtomId",   &System::maxAtomId)
            .def("maxBondId",   &System::maxBondId)
            .def("maxResidueId",&System::maxResidueId)
            .def("maxChainId",  &System::maxChainId)
            .def("maxCtId",     &System::maxCtId)

            /* list of element ids */
            .def("atoms",       list_atoms)
            .def("atomsAsList", list_atoms)
            .def("bonds",       list_bonds)
            .def("residues",    list_residues)
            .def("chains",      list_chains)
            .def("cts",         list_cts)


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
            .def("atomsForResidue", residue_atoms)
            .def("bondsForAtom",    atom_bonds)
            .def("residuesForChain",chain_residues)
            .def("chainsForCt",     ct_chains)
            .def("atomsForCt",      ct_atoms)
            .def("bondsForCt",      ct_bonds)

            /* tables */
            .def("tableNames",  table_names)
            .def("tableName",   &System::tableName)
            .def("table",       &System::table)
            .def("addTable",    wrap_system_add_table)
            .def("delTable",    &System::delTable)
            .def("removeTable", &System::removeTable)
            .def("renameTable", &System::renameTable)

            /* atom props */
            .def("setResidue",  &System::setResidue)
            .def("bondedAtoms", bonded_atoms)

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
            .def("selectAsList", list_Atomselect,
                    (arg("mol"),
                     arg("sel"),
                     arg("pos")=object(),
                     arg("box")=object()))
            .def("selectAsArray", array_Atomselect,
                    (arg("mol"),
                     arg("sel"),
                     arg("pos")=object(),
                     arg("box")=object()))

            /* append */
            .def("append", append_system)
            .def("clone",  wrap_clone)

            /* miscellaneous */
            .def("orderedIds",    ordered_ids)
            .def("updateFragids", update_fragids)
            .def("findBond",    &System::findBond)
            .def("provenance",      sys_provenance)
            .def("coalesceTables",    &System::coalesceTables)
            .def("translate",       sys_translate)
            .def("findContactIds",  sys_find_contact_ids)
            .def("topology",        sys_topology,
                    return_value_policy<manage_new_object>())
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

            /* PyCapsule conversion */
            .def("asCapsule", python::system_as_capsule)
            .staticmethod("asCapsule")
            .def("fromCapsule", python::system_from_capsule)
            .staticmethod("fromCapsule")
            ;
    def("HashSystem", HashSystem);

    class_<HydrogenBond>("HydrogenBond", no_init)
        .def("__init__", make_constructor(
                      init_hbond
                    , default_call_policies()
                    , (arg("d"),
                       arg("a"),
                       arg("h")=object(),
                       arg("c")=object(),
                       arg("ca")=object())),
                "Constructor: supply donor, acceptor and hydrogen positions")
        .add_property("energy", &HydrogenBond::energy,
                "Stride hbond energy function")
        .def_readwrite("r", &HydrogenBond::r,
                "donor-acceptor distance")
        .def_readwrite("p", &HydrogenBond::p,
                "donor-hydrogen-acceptor angle in radians")
        .def_readwrite("ti",&HydrogenBond::ti,
                "out-of-plane deviation of H from lone-pair plane")
        .def_readwrite("to",&HydrogenBond::to,
                "within-plane deviation of H from lone-pair bisector")
        ;
    }


}}
