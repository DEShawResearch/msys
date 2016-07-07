#include "io.hxx"
#include "dms.hxx"
#include "mae.hxx"
#include "pdb.hxx"
#include "psf.hxx"
#include "amber.hxx"
#include "mol2.hxx"
#include "xyz.hxx"
#include "sdf.hxx"

#include <boost/algorithm/string/predicate.hpp>
#include <boost/algorithm/string.hpp>

using namespace desres::msys;

namespace {
    // true if path ends with any of the endings in expr, separated by ','
    bool match(std::string const& path, std::string const& expr) {
        std::vector<std::string> endings;
        boost::split(endings, expr, boost::is_any_of(","));
        for (auto const& ending : endings) {
            if (boost::ends_with(path, ending)) return true;
        }
        return false;
    }

    bool match_web(std::string const& path) {
        return path.size()==4 && 
               isdigit(path[0]) &&
               isalpha(path[1]) &&
               isalnum(path[2]) &&
               isalnum(path[3]);
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
        "WEBPDB",
        "PSF"
    };

    class DefaultIterator : public LoadIterator {
        SystemPtr mol;
    public:
        DefaultIterator(std::string const& path, FileFormat format,
                        bool structure_only, bool without_tables) {
            mol = LoadWithFormat(path, format, structure_only, without_tables);
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
        to_lower(path);
        const char* DMS = "dms,dms.gz";
        const char* MAE = "mae,mae.gz,maegz,maeff,maeff.gz,cms,cms.gz";
        const char* PDB = "pdb";
        const char* PRM = "prmtop,prm7";
        const char* MOL2= "mol2";
        const char* XYZ = "xyz";
        const char* SDF = "sdf,sdf.gz,sdfgz";
        const char* PSF = "psf";

        if (match(path, DMS)) return DmsFileFormat;
        if (match(path, MAE)) return MaeFileFormat;
        if (match(path, PDB)) return PdbFileFormat;
        if (match(path, PRM)) return ParmTopFileFormat;
        if (match(path,MOL2)) return Mol2FileFormat;
        if (match(path, XYZ)) return XyzFileFormat;
        if (match(path, SDF)) return SdfFileFormat;
        if (match(path, PSF)) return PsfFileFormat;
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
                             bool structure_only, bool without_tables) {
        SystemPtr m;
        switch (format) {
            case DmsFileFormat: 
                m=ImportDMS(path, structure_only, without_tables); 
                break;
            case MaeFileFormat: 
                m=ImportMAE(path, false, structure_only, without_tables); 
                break;
            case PdbFileFormat: 
                m=ImportPDB(path); 
                break;
            case ParmTopFileFormat: 
                m=ImportPrmTop(path, structure_only, without_tables); 
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
            case PsfFileFormat:
                m=ImportPSF(path);
                break;
            case WebPdbFileFormat:
                m=ImportWebPDB(path);
            default:
                ;
        }
        return m;
    }

    SystemPtr Load(std::string const& path, FileFormat* opt_format,
                                            bool structure_only,
                                            bool without_tables) {
        FileFormat format = GuessFileFormat(path);
        if (opt_format) *opt_format = format;
        return LoadWithFormat(path, format, structure_only, without_tables);
    }

    LoadIteratorPtr LoadIterator::create(std::string const& path,
                                         FileFormat* opt_format,
                                         bool structure_only,
                                         bool without_tables) {

        FileFormat format = opt_format ? *opt_format : UnrecognizedFileFormat;
        if (!format) format=GuessFileFormat(path);
        if (opt_format) *opt_format = format;
        switch(format) {
            default:
                return LoadIteratorPtr(
                        new DefaultIterator(path, format, structure_only,
                                                          without_tables));
            case Mol2FileFormat:
                    return Mol2Iterator(path);
            case MaeFileFormat:
                    return MaeIterator(path, structure_only);
            case SdfFileFormat:
                    return SdfIterator(path);
            case PdbFileFormat:
                    return PDBIterator(path);
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
                ExportPDB(mol, path,
                      ( flags & SaveOptions::Append ? PDBExport::Append : 0)
                    );
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
                MSYS_FAIL("No support for saving file '" << path << "' of type "
                        << FileFormatAsString(format));
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
