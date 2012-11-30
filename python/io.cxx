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
        ExportSDF(mol,ss);
        return ss.str();
    }
}

namespace desres { namespace msys { 

    void export_io() {
        
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
        def("ImportMOL2", ImportMol2);
        def("ImportMOL2Many", import_mol2_many);
        def("ExportMOL2", ExportMol2);
        def("ImportXYZ", ImportXYZ);
        def("ExportSDFBytes", export_sdf_bytes);
        def("Load", Load,
                (arg("path"),
                 arg("opt_format")=object()));

    }
}}
