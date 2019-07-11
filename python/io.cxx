#include "wrap_obj.hxx"
#include "schema.hxx"
#include "mae.hxx"
#include "dms.hxx"
#include "pdb.hxx"
#ifndef _MSC_VER
#include "sdf.hxx"
#include "smiles.hxx"
#endif
#include "mol2.hxx"
#include "xyz.hxx"
#include "io.hxx"
#include "amber.hxx"
#include "json.hxx"

namespace {

#if PY_MAJOR_VERSION >= 3
    auto py_as_bytes = PyBytes_FromStringAndSize;
#else
    auto py_as_bytes = PyString_FromStringAndSize;
#endif

    std::string format_json(SystemPtr mol) {
        return FormatJson(mol, Provenance());
    }

    PyObject* format_dms(SystemPtr mol) {
        std::string contents = FormatDMS(mol, Provenance());
        return py_as_bytes(contents.data(), contents.size());
    }

    SystemPtr import_mae_from_buffer(PyObject* obj, bool ignore_unrecognized,
                                                    bool structure_only) {
        Py_buffer view[1];
        if (PyObject_GetBuffer(obj, view, PyBUF_ND)) {
            throw_error_already_set();
        }
        std::shared_ptr<Py_buffer> ptr(view, PyBuffer_Release);
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
        std::shared_ptr<Py_buffer> ptr(view, PyBuffer_Release);
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

    LoadIteratorPtr load_iterator_create(std::string const& path,
                                         bool structure_only) {
        return LoadIterator::create(path, structure_only);
    }
    SystemPtr load_iterator_next(LoadIterator& iter) {
        return iter.next();
    }

    void save(SystemPtr mol, std::string const& path, Provenance const& prov,
              bool append, bool structure_only) {
        unsigned flags = 0;
        if (append) flags |= SaveOptions::Append;
        if (structure_only) flags |= SaveOptions::StructureOnly;
        Save(mol, path, prov, flags);
    }

    SystemPtr import_dms(const std::string& path, bool structure_only=false) {
        return ImportDMS(path, structure_only);
    }
    SystemPtr import_mae(const std::string& path, bool ignore_unrecognized,
                         bool structure_only) {
        return ImportMAE(path, ignore_unrecognized, structure_only);
    }
    SystemPtr import_prmtop(const std::string& path, bool structure_only) {
        return ImportPrmTop(path, structure_only);
    }

    SystemPtr load(std::string const& path, bool structure_only, 
                                            bool without_tables) {
        return Load(path, NULL, structure_only, without_tables);
    }

    std::string indexed_file_path(IndexedFileLoader const& L) {
        return L.path();
    }
    size_t indexed_file_size(IndexedFileLoader const& L) {
        return L.size();
    }
    SystemPtr indexed_file_at(IndexedFileLoader const& L, size_t i) {
        return L.at(i);
    }

    void export_mol2(SystemPtr mol, std::string const& path,
                     Provenance const& prov, object _ids,
                     unsigned flags) {
        IdList ids;
        for (unsigned i=0, n=len(_ids); i<n; i++) {
            ids.push_back(extract<Id>(_ids[i]));
        }
        ExportMol2(mol, path, prov, ids, flags);
    }

    list import_pdb_unit_cell(double A, double B, double C,
                              double alpha, double beta, double gamma) {
        double cell[9];
        ImportPDBUnitCell(A,B,C,alpha,beta,gamma,cell);
        list L;
        for (int i=0; i<3; i++) {
            list sub;
            for (int j=0; j<3; j++) {
                sub.append(cell[3*i+j]);
            }
            L.append(sub);
        }
        return L;
    }

}

namespace desres { namespace msys { 

    void export_io() {
        
        enum_<DMSExport::Flags>("DMSExportFlags")
            .value("Default",           DMSExport::Default)
            .value("Append",            DMSExport::Append)
            .value("StructureOnly",     DMSExport::StructureOnly)
            .value("Unbuffered",        DMSExport::Unbuffered)
            ;

        enum_<MaeExport::Flags>("MaeExportFlags")
            .value("Default",           MaeExport::Default)
            .value("StructureOnly",     MaeExport::StructureOnly)
            //.value("CompressForcefield",MaeExport::CompressForcefield)
            .value("Append",            MaeExport::Append)
            ;

        enum_<Mol2Export::Flags>("Mol2ExportFlags")
            .value("Default",           Mol2Export::Default)
            .value("Append",            Mol2Export::Append)
            .value("MOE",               Mol2Export::MOE)
            ;

        enum_<PDBExport::Flags>("PDBExportFlags")
            .value("Default",           PDBExport::Default)
            .value("Append",            PDBExport::Append)
            ;

        class_<LoadIterator, LoadIteratorPtr, boost::noncopyable>("LoadIterator", no_init)
            .def("create", load_iterator_create).staticmethod("create")
            .def("next", load_iterator_next)
            ;

        class_<IndexedFileLoader, std::shared_ptr<IndexedFileLoader>,
            boost::noncopyable>("IndexedFileLoader", no_init)
            .def("create", IndexedFileLoader::create)
            .staticmethod("create")
            .def("path", indexed_file_path)
            .def("size", indexed_file_size)
            .def("at",   indexed_file_at)
            ;

        def("ImportDMS", import_dms);
        def("ImportDMSFromBuffer", import_dms_from_buffer);
        def("ExportDMS", ExportDMS);
        def("FormatDMS", format_dms);
        def("ImportMAE", import_mae);
        def("ImportMAEFromBuffer", import_mae_from_buffer);
        def("ExportMAE", ExportMAE);
        def("ExportMAEContents", ExportMAEContents);
        def("ImportPDB", ImportPDB);
        def("ExportPDB", ExportPDB);
        def("FetchPDB",  FetchPDB);
        def("ImportPrmTop", import_prmtop);
        def("ImportCrdCoordinates", ImportCrdCoordinates);
        def("ImportPDBCoordinates", ImportPDBCoordinates);
        def("ImportPDBUnitCell", import_pdb_unit_cell);
        def("ImportMOL2", ImportMol2);
        def("ImportMOL2Many", import_mol2_many);
        def("ExportMOL2", export_mol2);
        def("ImportXYZ", ImportXYZ);
        def("Load", load);
        def("Save", save);
#ifndef _MSC_VER
        def("FromSmilesString", FromSmilesString);
        def("ParseSDF", SdfTextIterator);
        def("FormatSDF", FormatSdf);
#endif
        def("FormatJson", format_json);
        def("ParseJson", ParseJson);
    }
}}
