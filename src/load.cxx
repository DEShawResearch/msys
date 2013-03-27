#include "load.hxx"
#include "atomsel/regex.hxx"
#include "dms.hxx"
#include "mae.hxx"
#include "pdb.hxx"
#include "amber.hxx"
#include "mol2.hxx"
#include "xyz.hxx"
#include "sdf.hxx"

#include <boost/algorithm/string.hpp>
#include <pcre.h>

using namespace desres::msys;

namespace {
    bool match(std::string path, std::string const& expr) {
        std::vector<std::string> endings;
        boost::split(endings, expr, boost::is_any_of(","));
        boost::to_lower(path);
        for (unsigned i=0; i<endings.size(); i++) {
            std::string p(".+\\.");
            p += endings[i];
            if (regex_match(path,atomsel::Regex(p))) {
                return true;
            }
        }
        return false;
    }

    const char *format_names[] = {
        "UNRECOGNIZED",
        "DMS",
        "MAE",
        "PDB",
        "PARM7",
        "MOL2",
        "XYZ",
        "SDF"
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

    FileFormat GuessFileFormat(std::string const& path) {

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

}}
