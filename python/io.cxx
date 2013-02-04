#include "wrap_obj.hxx"
#include "schema.hxx"
#include "mae.hxx"
#include "dms.hxx"
#include "pdb.hxx"
#include "mol2.hxx"
#include "xyz.hxx"
#include "load.hxx"
#include "amber.hxx"
#include "sdf.hxx"

namespace {

    SystemPtr import_mae_from_buffer(PyObject* obj, bool ignore_unrecognized,
                                                    bool structure_only) {
        Py_buffer view[1];
        if (PyObject_GetBuffer(obj, view, PyBUF_ND)) {
            throw_error_already_set();
        }
        boost::shared_ptr<Py_buffer> ptr(view, PyBuffer_Release);
        const char* bytes = reinterpret_cast<const char *>(view->buf);
        return ImportMAEFromBytes(bytes, view->len, 
                                  ignore_unrecognized,
                                  structure_only);
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

    list import_mol2_many(std::string const& path) {
        std::vector<SystemPtr> mols = ImportMol2Many(path);
        list L;
        for (unsigned i=0; i<mols.size(); i++) {
            L.append(object(mols[i]));
        }
        return L;
    }

    std::string export_sdf_bytes(SystemPtr mol) {
        std::stringstream ss;
        ExportSdf(mol,ss);
        return ss.str();
    }

    LoadIteratorPtr load_iterator_create(std::string const& path,
                                         bool structure_only) {
        return LoadIterator::create(path, structure_only);
    }
    SystemPtr load_iterator_next(LoadIterator& iter) {
        return iter.next();
    }

    void export_mae_many(object ctarr, std::string const& path,
                         Provenance const& prov,
                         bool with_forcefield,
                         bool with_compression) {
        std::vector<SystemPtr> ctlist;
        for (Py_ssize_t i=0; i<len(ctarr); i++) {
            ctlist.push_back(extract<SystemPtr>(ctarr[i]));
        }
        ExportMAEMany(ctlist, path, prov, with_forcefield, with_compression);
    }

    void export_dms_many(object ctarr, std::string const& path,
                         Provenance const& prov) {
        std::vector<SystemPtr> ctlist;
        for (Py_ssize_t i=0; i<len(ctarr); i++) {
            ctlist.push_back(extract<SystemPtr>(ctarr[i]));
        }
        ExportDMSMany(ctlist, path, prov);
    }
}

namespace desres { namespace msys { 

    void export_io() {
        
        class_<LoadIterator, LoadIteratorPtr, boost::noncopyable>("LoadIterator", no_init)
            .def("create", load_iterator_create).staticmethod("create")
            .def("next", load_iterator_next)
            ;

        def("ImportDMS", ImportDMS);
        def("ImportDMSFromBuffer", import_dms_from_buffer);
        def("ExportDMS", ExportDMS);
        def("ExportDMSMany", export_dms_many);
        def("ImportMAE", ImportMAE);
        def("ImportMAEFromBuffer", import_mae_from_buffer);
        def("ExportMAE", ExportMAE);
        def("ExportMAEMany", export_mae_many);
        def("ImportPDB", ImportPDB);
        def("ExportPDB", ExportPDB);
        def("ImportPrmTop", ImportPrmTop);
        def("ImportCrdCoordinates", ImportCrdCoordinates);
        def("ImportMOL2", ImportMol2);
        def("ImportMOL2Many", import_mol2_many);
        def("ExportMOL2", ExportMol2);
        def("ImportXYZ", ImportXYZ);
        def("ExportSDFBytes", export_sdf_bytes);
        def("Load", Load,
                (arg("path"),
                 arg("structure_only")=false,
                 arg("opt_format")=object()));

    }
}}
