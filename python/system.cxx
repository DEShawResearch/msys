#include "pymod.hxx"
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <numpy/ndarrayobject.h>
#include "capsule.hxx"

#include <msys/import.hxx>
#include <msys/schema.hxx>
#include <msys/atomsel.hxx>
#include <msys/append.hxx>
#include <msys/clone.hxx>
#include <msys/geom.hxx>
#include <msys/dms.hxx>
#include <msys/hbond.hxx>
#include <msys/contacts.hxx>
#include <msys/hash.hxx>

#include <pfx/graph.hxx>
#include <pfx/cell.hxx>

namespace {

    handle atom_prop_type(System& sys, Id col) {
        return from_value_type(sys.atomPropType(col));
    }
    Id add_atom_prop( System& sys, String const& name, object type ) {
        return sys.addAtomProp(name, as_value_type(type));
    }

    handle bond_prop_type(System& sys, Id col) {
        return from_value_type(sys.bondPropType(col));
    }
    Id add_bond_prop( System& sys, String const& name, object type ) {
        return sys.addBondProp(name, as_value_type(type));
    }

    MultiIdList update_fragids(System& p) {
        MultiIdList fragments;
        p.updateFragids(&fragments);
        return fragments;
    }

    Provenance prov_from_args(std::vector<std::string> args) {
        if (args.empty()) return Provenance();
        std::vector<char *> argv;
        for (auto& s : args) argv.push_back(const_cast<char *>(s.data()));
        return Provenance::fromArgs(argv.size(), argv.data());
    }

    std::vector<Provenance> sys_provenance( System const& sys ) {
        return sys.provenance();
    }
    void sys_set_provenance(System& sys, std::vector<Provenance>& prov) {
        sys.provenance().swap(prov);
    }

    array global_cell(object sysobj) {
        auto sys = sysobj.cast<System*>();
        return array_t<double>({3,3}, sys->global_cell[0], sysobj);
    }

    void getpos3(object x, double *a) {
        auto arr = PyArray_FromAny(
                x.ptr(),
                PyArray_DescrFromType(NPY_DOUBLE),
                1, 1,
                NPY_ARRAY_C_CONTIGUOUS | NPY_FORCECAST,
                nullptr);
        if (!arr) throw error_already_set();
        memcpy(a, PyArray_DATA(arr), 3*sizeof(double));
        Py_DECREF(arr);
    }

    double py_calc_distance(object A, object B) {
        double a[3], b[3];
        getpos3(A, a);
        getpos3(B, b);
        return calc_distance(a,b);
    }
    double py_calc_vec_angle(object A, object B) {
        double a[3], b[3];
        getpos3(A, a);
        getpos3(B, b);
        return calc_vec_angle(a,b);
    }
    double py_calc_angle(object A, object B, object C) {
        double a[3], b[3], c[3];
        getpos3(A, a);
        getpos3(B, b);
        getpos3(C, c);
        return calc_angle(a,b,c);
    }
    double py_calc_vec_dihedral(object A, object B, object C) {
        double a[3], b[3], c[3];
        getpos3(A, a);
        getpos3(B, b);
        getpos3(C, c);
        return calc_vec_dihedral(a,b,c);
    }
    double py_calc_dihedral(object A, object B, object C, object D) {
        double a[3], b[3], c[3], d[3];
        getpos3(A, a);
        getpos3(B, b);
        getpos3(C, c);
        getpos3(D, d);
        return calc_dihedral(a,b,c,d);
    }

    array py_apply_dihedral_geometry(object A, object B, object C,
                                      double r, double theta, double phi) {
        double a[3], b[3], c[3];
        getpos3(A, a);
        getpos3(B, b);
        getpos3(C, c);
        auto arr = array_t<double>(3);
        apply_dihedral_geometry(arr.mutable_data(),a,b,c,r,theta,phi);
        return arr;
    }

    double py_calc_planarity(array_t<double, array::c_style | array::forcecast> arr) {
        if (arr.ndim()!=2 || arr.shape(1)!=3) {
            PyErr_Format(PyExc_ValueError, 
                    "Supplied %ld-d positions, expected 3-d", 
                    arr.shape(1));
            throw error_already_set();
        }
        return calc_planarity(arr.shape(0), arr.data());
    }

    bool py_line_intersects_tri(object A, object B, object C, 
                                object R, object S) {
        double a[3], b[3], c[3], r[3], s[3];
        getpos3(A, a);
        getpos3(B, b);
        getpos3(C, c);
        getpos3(R, r);
        getpos3(S, s);
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
        if (sys.maxAtomId() != sys.atomCount()) {
            PyErr_Format(PyExc_ValueError, "System has deleted atoms");
            throw error_already_set();
        }

        const Id* idptr=NULL, *otherptr=NULL;
        unsigned nids=0, nother=0;
        const double* posptr=NULL;
        std::vector<double> posdata;
        IdList ids;
        array_t<unsigned> idarr, otherarr;
        array_t<double> posarr;
        bool skip_symmetric = true;

        if (idobj.is_none()) {
            ids = sys.atoms();
            idptr = ids.data();
            nids = ids.size();
        } else {
            idarr = array_t<unsigned>::ensure(idobj);
            if (!idarr) throw error_already_set();
            idptr = idarr.data();
            nids = idarr.size();
        }

        if (otherobj.is_none()) {
            otherptr = idptr;
            nother=nids;
        } else {
            skip_symmetric = false;
            otherarr = array_t<unsigned>::ensure(otherobj);
            if (!otherarr) throw error_already_set();
            otherptr = otherarr.data();
            nother = otherarr.size();
        }

        if (posobj.is_none()) {
            sys.getPositions(std::back_inserter(posdata));
            posptr = posdata.data();
        } else {
            posarr = array_t<double>::ensure(posobj);
            posptr = posarr.data();
        }

        list result;
        find_contacts(cutoff, posptr,
                      idptr, idptr+nids,
                      otherptr, otherptr+nother,
                      Finder(sys, skip_symmetric, result));
        return result;
    }


    array sys_getpos(System const& sys, object idobj) {
        array_t<double> arr;
        if (idobj.is_none()) {
            size_t dims[] = {sys.atomCount(), 3};
            arr = array_t<double>(dims);
            sys.getPositions(arr.mutable_data());
        } else {
            auto idarr = array_t<unsigned>::ensure(idobj);
            if (!idarr) throw error_already_set();
            ssize_t dims[] = {idarr.size(), 3};
            arr = array_t<double>(dims);
            auto ids = idarr.data();
            auto ptr = arr.mutable_data();
            for (int i=0, n=idarr.size(); i<n; i++) {
                auto id = ids[i];
                if (!sys.hasAtom(id)) {
                    PyErr_Format(PyExc_ValueError, "Invalid id %u", id);
                    throw error_already_set();
                }
                atom_t const& atm = sys.atomFAST(id);
                ptr[0] = atm.x;
                ptr[1] = atm.y;
                ptr[2] = atm.z;
                ptr += 3;
            }
        }
        return arr;
    }

    array sys_getvel(System const& sys, object idobj) {
        array_t<double> arr;
        if (idobj.is_none()) {
            size_t dims[] = {sys.atomCount(), 3};
            arr = array_t<double>(dims);
            sys.getVelocities(arr.mutable_data());
        } else {
            auto idarr = array_t<unsigned>::ensure(idobj);
            if (!idarr) throw error_already_set();
            size_t dims[] = {sys.atomCount(), 3};
            arr = array_t<double>(dims);
            auto ids = idarr.data();
            auto ptr = arr.mutable_data();
            for (int i=0, n=idarr.size(); i<n; i++) {
                auto id = ids[i];
                if (!sys.hasAtom(id)) {
                    PyErr_Format(PyExc_ValueError, "Invalid id %u", id);
                    throw error_already_set();
                }
                atom_t const& atm = sys.atomFAST(id);
                ptr[0] = atm.vx;
                ptr[1] = atm.vy;
                ptr[2] = atm.vz;
                ptr += 3;
            }
        }
        return arr;
    }

    void sys_setpos(System& sys, array_t<double, array::c_style | array::forcecast> arr, object idobj) {
        if (arr.ndim()!=2 || arr.shape(1)!=3) {
            PyErr_Format(PyExc_ValueError, 
                    "Supplied %ld-d positions, expected 3-d", 
                    arr.shape(1));
            throw error_already_set();
        }
        auto pos = arr.data();
        auto n = arr.shape(0);

        if (idobj.is_none()) {
            if (n != sys.atomCount()) {
                PyErr_Format(PyExc_ValueError, 
                        "Supplied %u positions, but system has %u atoms", 
                        n, sys.atomCount());
                throw error_already_set();
            }
            sys.setPositions(pos);
        } else {
            auto idarr = array_t<unsigned>::ensure(idobj);
            if (!idarr) throw error_already_set();
            if (n != idarr.shape(0)) {
                PyErr_Format(PyExc_ValueError,
                        "Supplied %u positions != %u ids",
                        n, idarr.shape(0));
                throw error_already_set();
            }
            auto ids = idarr.data();
            for (int i=0; i<n; i++) {
                Id id = ids[i];
                if (!sys.hasAtom(id)) {
                    PyErr_Format(PyExc_ValueError,
                            "Id %u at position %ld is invalid", id, i);
                    throw error_already_set();
                }
                atom_t& atm = sys.atomFAST(id);
                atm.x = pos[0];
                atm.y = pos[1];
                atm.z = pos[2];
                pos += 3;
            }
        }
    }

    void sys_setvel(System& sys, array_t<double, array::c_style | array::forcecast> arr, object idobj) {
        if (arr.ndim()!=2 || arr.shape(1)!=3) {
            PyErr_Format(PyExc_ValueError, 
                    "Supplied %ld-d velocities, expected 3-d", 
                    arr.shape(1));
            throw error_already_set();
        }
        auto pos = arr.data();
        auto n = arr.shape(0);

        if (idobj.is_none()) {
            if (n != sys.atomCount()) {
                PyErr_Format(PyExc_ValueError, 
                        "Supplied %u velocities, but system has %u atoms", 
                        n, sys.atomCount());
                throw error_already_set();
            }
            sys.setVelocities(pos);
        } else {
            auto idarr = array_t<unsigned>::ensure(idobj);
            if (!idarr) throw error_already_set();
            if (n != idarr.shape(0)) {
                PyErr_Format(PyExc_ValueError,
                        "Supplied %u velocities != %u ids",
                        n, idarr.shape(0));
                throw error_already_set();
            }
            auto ids = idarr.data();
            for (int i=0; i<n; i++) {
                Id id = ids[i];
                if (!sys.hasAtom(id)) {
                    PyErr_Format(PyExc_ValueError,
                            "Id %u at velocity %ld is invalid", id, i);
                    throw error_already_set();
                }
                atom_t& atm = sys.atomFAST(id);
                atm.vx = pos[0];
                atm.vy = pos[1];
                atm.vz = pos[2];
                pos += 3;
            }
        }
    }

    IdList wrap_atomselect(SystemPtr mol, std::string const& sel,
                           object posobj, object boxobj) {
        array_t<float> posarr;
        array_t<double> boxarr;
        const float* pos = NULL;
        const double* box = NULL;
        if (!posobj.is_none()) {
            posarr = array_t<float>::ensure(posobj);
            if (!posarr) throw error_already_set();
            if (posarr.ndim()!=2 || posarr.shape(0)!=mol->atomCount() || posarr.shape(1)!=3) {
                PyErr_Format(PyExc_ValueError, "pos has wrong shape");
                throw error_already_set();
            }
            pos = posarr.data();
        }
        if (!boxobj.is_none()) {
            boxarr = array_t<double>::ensure(boxobj);
            if (!boxarr) throw error_already_set();
            if (boxarr.ndim()!=2 || boxarr.shape(0)!=3 || boxarr.shape(1)!=3) {
                PyErr_Format(PyExc_ValueError, "box has wrong shape");
                throw error_already_set();
            }
            box = boxarr.data();
        }
        return Atomselect(mol, sel, pos, box);
    }

    object array_Atomselect(SystemPtr mol, std::string const& sel,
                               object pos, object box) {
        auto ids = wrap_atomselect(mol,sel, pos, box);
        auto arr = array_t<unsigned>(ids.size());
        memcpy(arr.mutable_data(), ids.data(), ids.size()*sizeof(ids[0]));
        return arr;
    }

    const double* get_vec3d(object obj) {
        if (obj.is_none()) return nullptr;
        auto arr = array_t<double>::ensure(obj);
        if (!arr) return nullptr;
        if (arr.size()!=3) {
            PyErr_Format(PyExc_ValueError,
                    "Expected 3 elements in vector, got %ld",
                    arr.size());
            throw error_already_set();
        }
        return arr.data();
    }

    HydrogenBond *init_hbond(
            object dobj,
            object aobj,
            object hobj,
            object cobj,
            object caobj) {
        return new HydrogenBond(
                get_vec3d(dobj),
                get_vec3d(aobj),
                get_vec3d(hobj),
                get_vec3d(cobj),
                get_vec3d(caobj));
    }

    pfx::Graph* sys_topology(SystemPtr mol) { return new pfx::Graph(mol); }

    TermTablePtr wrap_system_add_table(SystemPtr mol, std::string const& name, Id natoms, object obj) {
        ParamTablePtr params;
        if (!obj.is_none()) params = obj.cast<ParamTablePtr>();
        return mol->addTable(name, natoms, params);
    }

    SystemPtr init_from_pickle(std::string const& format, buffer bobj) {
        if (format=="dmscontents") {
            auto r = bobj.request();
            if (r.size < 100) {
                throw std::runtime_error("pickle contents are too short: len = " + std::to_string(r.size));
            }
            auto buf = reinterpret_cast<unsigned char *>(r.ptr);
            // check for zlib-compressed data.  
            if (buf[0] == 0x78 && (
                        // alternatives based on possible zlib compression levels
                        buf[1] == 0x01 ||
                        buf[1] == 0x5e ||
                        buf[1] == 0x9c ||
                        buf[1] == 0xda)) {

                auto decompress = module::import("zlib").attr("decompress");
                auto dobj = decompress(bobj);
                auto dbuf = buffer(dobj);
                auto dreq = dbuf.request();
                return ImportDMSFromBytes((const char *)dreq.ptr, dreq.size);
            }
            return ImportDMSFromBytes((const char *)r.ptr, r.size);
        }
        throw std::runtime_error("Unsupported pickle format " + format);
    }

}

namespace desres { namespace msys { 

    void export_system(module m) {
        _import_array();
        if (PyErr_Occurred()) throw error_already_set();

        m.def("bad", bad);
        m.attr("BadId") = cast((Id)BadId);

        m.def("calc_distance", py_calc_distance);
        m.def("calc_vec_angle", py_calc_vec_angle);
        m.def("calc_angle", py_calc_angle);
        m.def("calc_vec_dihedral", py_calc_vec_dihedral);
        m.def("calc_dihedral", py_calc_dihedral);
        m.def("calc_planarity", py_calc_planarity);
        m.def("line_intersects_tri", py_line_intersects_tri);
        m.def("apply_dihedral_geometry", py_apply_dihedral_geometry);

        class_<pfx::Graph>(m, "Topology")
            .def(init<unsigned>())
            .def_property_readonly("nverts", &pfx::Graph::nverts)
            .def_property_readonly("nedges", &pfx::Graph::nedges)
            ;

        class_<NonbondedInfo>(m, "NonbondedInfo")
            .def_readwrite("vdw_funct", &NonbondedInfo::vdw_funct,
                    "Name of the vdw functional form; e.g., 'vdw_12_6'")
            .def_readwrite("vdw_rule", &NonbondedInfo::vdw_rule,
                    "Nonbonded combining rule; e.g., 'arithmetic/geometric'")
            .def_readwrite("es_funct", &NonbondedInfo::es_funct,
                    "Name of the electrostatic functional form")
            ;

        class_<Provenance>(m, "Provenance")
            .def(init<>())
            .def_static("fromArgs", prov_from_args)
            .def_readwrite("version", &Provenance::version)
            .def_readwrite("timestamp", &Provenance::timestamp)
            .def_readwrite("user", &Provenance::user)
            .def_readwrite("workdir", &Provenance::workdir)
            .def_readwrite("cmdline", &Provenance::cmdline)
            .def_readwrite("executable", &Provenance::executable)
            ;

        enum_<CloneOption::Flags>(m, "CloneOption", arithmetic())
            .value("Default",       CloneOption::Default)
            .value("ShareParams",   CloneOption::ShareParams)
            .value("UseIndex",      CloneOption::UseIndex)
            .value("StructureOnly", CloneOption::StructureOnly)
            ;

        m.def("TableSchemas", TableSchemas);
        m.def("NonbondedSchemas", NonbondedSchemas);

        class_<System,SystemPtr>(m, "SystemPtr")
            .def("__eq__", [](System* self, System* other) { return self==other; })
            .def("__ne__", [](System* self, System* other) { return self!=other; })
            .def("__hash__", [](System* g) { return size_t(g); })
            .def_static("create",  &System::create)
            .def_readwrite("name", &System::name)
            .def_property_readonly("global_cell", global_cell)
            .def_readwrite("nonbonded_info", &System::nonbonded_info)

            /* pickle support */
            .def(pickle(
                [](SystemPtr mol) { // __getstate__
                    auto bobj = bytes(FormatDMS(mol, Provenance()));
                    auto compress = module::import("zlib").attr("compress");
                    auto zobj = compress(bobj, 1);
                    return make_tuple("dmscontents", zobj);
                },
                [](tuple t) { // __setstate__
                    if (t.size() != 2) throw std::runtime_error("Invalid state!");
                    return init_from_pickle(t[0].str(), object(t[1]));
                    }
                ))

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
            .def("atoms",       &System::atoms)
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
            .def("atomsForResidue", &System::atomsForResidue)
            .def("bondsForAtom",    &System::bondsForAtom)
            .def("residuesForChain",&System::residuesForChain)
            .def("chainsForCt",     &System::chainsForCt)
            .def("atomsForCt",      &System::atomsForCt)
            .def("bondsForCt",      &System::bondsForCt)

            /* tables */
            .def("tableNames",  &System::tableNames)
            .def("tableName",   &System::tableName)
            .def("table",       &System::table)
            .def("addTable",    wrap_system_add_table)
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

            /* extended bond props */
            .def("bondPropCount",&System::bondPropCount)
            .def("bondPropName", &System::bondPropName)
            .def("bondPropIndex",&System::bondPropIndex)
            .def("bondPropType", bond_prop_type)
            .def("addBondProp",  add_bond_prop)
            .def("delBondProp",  &System::delBondProp)

            /* auxiliary tables */
            .def("auxTableNames",&System::auxTableNames)
            .def("auxTable",     &System::auxTable)
            .def("addAuxTable",  &System::addAuxTable)
            .def("delAuxTable",  &System::delAuxTable)
            .def("removeAuxTable",&System::removeAuxTable)

            /* schemas */
            .def("addTableFromSchema", AddTable)
            .def("addNonbondedFromSchema", AddNonbonded)

            /* atom selection */
            .def("selectAsList", wrap_atomselect, arg("sel"), arg("pos")=none(), arg("box")=none())
            .def("selectAsArray", array_Atomselect, arg("sel"), arg("pos")=none(), arg("box")=none())

            /* append */
            .def("append", AppendSystem)
            .def("clone",  Clone)

            /* miscellaneous */
            .def("orderedIds",    &System::orderedIds)
            .def("updateFragids", update_fragids)
            .def("findBond",    &System::findBond)
            .def("provenance",      sys_provenance)
            .def("setProvenance", sys_set_provenance)
            .def("coalesceTables",    &System::coalesceTables)
            .def("translate",       sys_translate)
            .def("findContactIds",  sys_find_contact_ids)
            .def("topology",        sys_topology)
            .def("getPositions", sys_getpos, arg("ids")=none())
            .def("setPositions",    sys_setpos, arg("pos"), arg("ids")=none())
            .def("getVelocities", sys_getvel, arg("ids")=none())
            .def("setVelocities",    sys_setvel, arg("vel"), arg("ids")=none())

            /* PyCapsule conversion */
            .def_static("asCapsule", [](SystemPtr ptr) -> handle { return python::system_as_capsule(ptr); })
            .def_static("fromCapsule", [](handle h) { return python::system_from_capsule(h.ptr()); })
            .def("addProvenance", &System::addProvenance)
            ;
    m.def("HashSystem", HashSystem);

    class_<SystemImporter>(m, "SystemImporter")
        .def(init<SystemPtr>())
        .def("initialize", &SystemImporter::initialize)
        .def("terminateChain", &SystemImporter::terminateChain)
        .def("addAtom", &SystemImporter::addAtom)
        ;


    class_<HydrogenBond>(m, "HydrogenBond", dynamic_attr())
        .def(init(&init_hbond), 
                    arg("d"),
                    arg("a"),
                    arg("h")=none(),
                    arg("c")=none(),
                    arg("ca")=none())
        .def_property_readonly("energy", &HydrogenBond::energy,
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
