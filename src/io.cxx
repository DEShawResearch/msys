#include "io.hxx"
#include "atomsel/regex.hxx"
#include "dms.hxx"
#include "mae.hxx"
#include "pdb.hxx"
#include "amber.hxx"
#include "mol2.hxx"
#include "xyz.hxx"
#include "sdf.hxx"

#include <boost/algorithm/string.hpp>

using namespace desres::msys;

namespace {
    bool match(std::string const& path, std::string const& expr) {
        std::vector<std::string> endings;
        boost::split(endings, expr, boost::is_any_of(","));
        for (unsigned i=0; i<endings.size(); i++) {
            std::string p(".+\\.");
            p += endings[i];
            if (regex_match(path,atomsel::Regex(p))) {
                return true;
            }
        }
        return false;
    }

    bool match_web(std::string const& path) {
        return regex_match(path, atomsel::Regex("[0-9][a-z][a-z][a-z]"));
    }

    const char *format_names[] = {
        "UNRECOGNIZED",
        "DMS",
        "MAE",
        "PDB",
        "PARM7",
        "MOL2",
        "XYZ",
        "SDF",
        "WEBPDB"
    };

    class DefaultIterator : public LoadIterator {
        SystemPtr mol;
    public:
        DefaultIterator(std::string const& path, FileFormat format,
                        bool structure_only) {
            mol = LoadWithFormat(path, format, structure_only);
        }
        SystemPtr next() {
            /* returns original mol the first time, then NULL */
            SystemPtr ptr;
            ptr.swap(mol);
            return ptr;
        }
    };
}

namespace desres { namespace msys {

    FileFormat GuessFileFormat(std::string const& _path) {

        std::string path(_path);
        boost::to_lower(path);
        const char* DMS = "dms,dms.gz";
        const char* MAE = "mae,mae.gz,maegz,maeff,maeff.gz,cms,cms.gz";
        const char* PDB = "pdb";
        const char* PRM = "prmtop,prm7";
        const char* MOL2= "mol2";
        const char* XYZ = "xyz";
        const char* SDF = "sdf,sdf.gz,sdfgz";

        if (match(path, DMS)) return DmsFileFormat;
        if (match(path, MAE)) return MaeFileFormat;
        if (match(path, PDB)) return PdbFileFormat;
        if (match(path, PRM)) return ParmTopFileFormat;
        if (match(path,MOL2)) return Mol2FileFormat;
        if (match(path, XYZ)) return XyzFileFormat;
        if (match(path, SDF)) return SdfFileFormat;
        if (match_web(path))  return WebPdbFileFormat;
        return UnrecognizedFileFormat;
    }

    std::string FileFormatAsString(FileFormat format) {
        return format_names[format];
    }

    FileFormat FileFormatFromString(std::string const& name) {
        unsigned i,n = sizeof(format_names)/sizeof(format_names[0]);
        for (i=0; i<n; i++) {
            if (name==format_names[i]) {
                return FileFormat(i);
            }
        }
        return UnrecognizedFileFormat;
    }

    SystemPtr LoadWithFormat(std::string const& path, FileFormat format,
                             bool structure_only) {
        SystemPtr m;
        switch (format) {
            case DmsFileFormat: 
                m=ImportDMS(path, structure_only); 
                break;
            case MaeFileFormat: 
                m=ImportMAE(path, false, structure_only); 
                break;
            case PdbFileFormat: 
                m=ImportPDB(path); 
                break;
            case ParmTopFileFormat: 
                m=ImportPrmTop(path, structure_only); 
                break;
            case Mol2FileFormat: 
                m=ImportMol2(path); 
                break;
            case XyzFileFormat: 
                m=ImportXYZ(path); 
                break;
            case SdfFileFormat:
                m=ImportSdf(path);
                break;
            case WebPdbFileFormat:
                m=ImportWebPDB(path);
            default:
                ;
        }
        return m;
    }

    SystemPtr Load(std::string const& path, bool structure_only,
                   FileFormat* opt_format) {
        FileFormat format = GuessFileFormat(path);
        if (opt_format) *opt_format = format;
        return LoadWithFormat(path, format, structure_only);
    }

    LoadIteratorPtr LoadIterator::create(std::string const& path,
                                         bool structure_only,
                                         FileFormat* opt_format) {

        FileFormat format = opt_format ? *opt_format : UnrecognizedFileFormat;
        if (!format) format=GuessFileFormat(path);
        if (opt_format) *opt_format = format;
        switch(format) {
            default:
                return LoadIteratorPtr(
                        new DefaultIterator(path, format, structure_only));
            case Mol2FileFormat:
                    return Mol2Iterator(path);
            case MaeFileFormat:
                    return MaeIterator(path, structure_only);
            case SdfFileFormat:
                    return SdfIterator(path);
            case UnrecognizedFileFormat:
                MSYS_FAIL("Unable to determine format of '" << path << "'");
        }
    }

    void SaveWithFormat(SystemPtr mol, 
                        std::string const& path, 
                        Provenance const& prov,
                        FileFormat format,
                        unsigned flags) {

        switch (format) {
            case DmsFileFormat: 
                ExportDMS(mol, path, prov, 
                    ( flags & SaveOptions::Append ? DMSExport::Append : 0)
                  | (flags & SaveOptions::StructureOnly ? DMSExport::StructureOnly : 0)
                    );
                break;
            case MaeFileFormat: 
                ExportMAE(mol, path, prov, 
                      (flags & SaveOptions::Append ? MaeExport::Append : 0) 
                    | (flags & SaveOptions::StructureOnly ? MaeExport::StructureOnly : 0)
                    );
                break;
            case PdbFileFormat: 
                if (flags != 0) {
                    MSYS_FAIL("PDB export supports only default save option");
                }
                ExportPDB(mol, path);
                break;
            case ParmTopFileFormat: 
                MSYS_FAIL("PRM/TOP export not supported");
                break;
            case Mol2FileFormat: 
                ExportMol2(mol, path, prov,
                      (flags & SaveOptions::Append ? Mol2Export::Append : 0)
                    );
                break;
            case XyzFileFormat: 
                if (flags != 0) {
                    MSYS_FAIL("XYZ export supports only default save option");
                }
                ExportXYZ(mol, path);
                break;
            case SdfFileFormat:
                ExportSdf(mol, path,
                      (flags & SaveOptions::Append ? SdfExport::Append : 0)
                    );
                break;
            default:
                ;
        }
    }

    void Save(SystemPtr mol, 
              std::string const& path, 
              Provenance const& prov,
              unsigned flags) {
        FileFormat fmt = GuessFileFormat(path);
        if (fmt == UnrecognizedFileFormat) {
            MSYS_FAIL("Unable to determine format of '" << path << "'");
        }
        SaveWithFormat(mol, path, prov, fmt, flags);
    }


}}
