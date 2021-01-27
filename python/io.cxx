#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

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

using namespace pybind11;
using namespace desres::msys;

namespace {

    std::string format_json(SystemPtr mol, int maxDecimals) {
        unsigned flags = 0;
        return FormatJson(mol, Provenance(), flags, maxDecimals);
    }

    bytes format_dms(SystemPtr mol) {
        return FormatDMS(mol, Provenance());
    }

    SystemPtr import_mae_from_buffer(buffer buf, bool ignore_unrecognized,
                                                 bool structure_only) {

        auto req = buf.request();
        return ImportMAEFromBytes((const char*)req.ptr, req.size,
                                  ignore_unrecognized,
                                  structure_only);
    }

    SystemPtr import_dms_from_buffer(buffer buf, bool structure_only ) {
        auto req = buf.request();
        return ImportDMSFromBytes((const char*)req.ptr, req.size, structure_only);
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


    list import_pdb_unit_cell(double A, double B, double C,
                              double alpha, double beta, double gamma) {
        double cell[9];
        ImportPDBUnitCell(A,B,C,alpha,beta,gamma,cell);
        list L;
        for (int i=0; i<3; i++) {
            list sub;
            for (int j=0; j<3; j++) {
                sub.append(cast(cell[3*i+j]));
            }
            L.append(sub);
        }
        return L;
    }

    object format_sdf(SystemPtr mol, bool as_bytes) {
        std::string s(FormatSdf(mol));
        if (as_bytes) {
            return bytes(s);
        }
        auto ptr = PyUnicode_FromStringAndSize(s.data(), s.size());
        if (!ptr) {
            PyErr_SetString(PyExc_ValueError, "system could not be UTF8-encoded; try calling with as_bytes=True");
            throw error_already_set();
        }
        return reinterpret_steal<str>(ptr);
    }
}

namespace desres { namespace msys { 

    void export_io(module m) {
        
        enum_<DMSExport::Flags>(m, "DMSExportFlags", arithmetic())
            .value("Default",           DMSExport::Default)
            .value("Append",            DMSExport::Append)
            .value("StructureOnly",     DMSExport::StructureOnly)
            .value("Unbuffered",        DMSExport::Unbuffered)
            ;

        enum_<MaeExport::Flags>(m, "MaeExportFlags", arithmetic())
            .value("Default",           MaeExport::Default)
            .value("StructureOnly",     MaeExport::StructureOnly)
            //.value("CompressForcefield",MaeExport::CompressForcefield)
            .value("Append",            MaeExport::Append)
            .value("AllowReorderAtoms", MaeExport::AllowReorderAtoms)
            ;

        enum_<Mol2Export::Flags>(m, "Mol2ExportFlags", arithmetic())
            .value("Default",           Mol2Export::Default)
            .value("Append",            Mol2Export::Append)
            .value("MOE",               Mol2Export::MOE)
            ;

        enum_<PDBExport::Flags>(m, "PDBExportFlags", arithmetic())
            .value("Default",           PDBExport::Default)
            .value("Append",            PDBExport::Append)
            .value("Reorder",           PDBExport::Reorder)
            ;

        class_<LoadIterator, LoadIteratorPtr>(m, "LoadIterator")
            .def_static("create", [](std::string const& path, bool structure_only) { return LoadIterator::create(path, structure_only); })
            .def("next", [](LoadIterator& li) { return li.next(); })
            ;

        class_<IndexedFileLoader, std::shared_ptr<IndexedFileLoader>>(m, "IndexedFileLoader")
            .def_static("create", &IndexedFileLoader::create)
            .def("path", [](IndexedFileLoader& self) { return self.path(); })
            .def("size", [](IndexedFileLoader& self) { return self.size(); })
            .def("at", [](IndexedFileLoader& self, size_t entry) { return self.at(entry); })
            ;

        m.def("ImportDMS", import_dms);
        m.def("ImportDMSFromBuffer", import_dms_from_buffer);
        m.def("ExportDMS", ExportDMS);
        m.def("FormatDMS", format_dms);
        m.def("ImportMAE", import_mae);
        m.def("ImportMAEFromBuffer", import_mae_from_buffer);
        m.def("ExportMAE", ExportMAE);
        m.def("ExportMAEContents", ExportMAEContents);
        m.def("ImportPDB", ImportPDB);
        m.def("ExportPDB", ExportPDB);
        m.def("FetchPDB",  FetchPDB);
        m.def("ImportPrmTop", import_prmtop);
        m.def("ImportCrdCoordinates", ImportCrdCoordinates);
        m.def("ImportPDBCoordinates", ImportPDBCoordinates);
        m.def("ImportPDBUnitCell", import_pdb_unit_cell);
        m.def("ImportMOL2", ImportMol2);
        m.def("ImportMOL2Many", ImportMol2Many);
        m.def("ExportMOL2", ExportMol2);
        m.def("ImportXYZ", ImportXYZ);
        m.def("Load", load);
        m.def("Save", save);
        m.def("FromSmilesString", FromSmilesString);
        m.def("ParseSDF", SdfTextIterator);
        m.def("FormatSDF", format_sdf, arg("system"), arg("as_bytes")=false);
        m.def("FormatJson", format_json);
        m.def("ParseJson", ParseJson);
    }
}}
