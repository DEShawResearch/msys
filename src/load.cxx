#include "load.hxx"
#include "atomsel/regex.hxx"
#include "dms.hxx"
#include "mae.hxx"
#include "pdb.hxx"
#include "amber.hxx"

#include <boost/algorithm/string.hpp>
#include <pcre.h>

using namespace desres::msys;

namespace {
    bool match(std::string const& path, std::string const& expr) {
        std::vector<std::string> endings;
        boost::split(endings, expr, boost::is_any_of(","));
        for (unsigned i=0; i<endings.size(); i++) {
            std::string p(".+\\.");
            p += endings[i];
            atomsel::Regex regex(p, PCRE_CASELESS);
            if (regex.match(path)) {
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
        "PARM7"
    };
}

namespace desres { namespace msys {

    FileFormat GuessFileFormat(std::string const& path) {

        const char* DMS = "dms,dms.gz";
        const char* MAE = "mae,mae.gz,maegz,maeff,maeff.gz,cms,cms.gz";
        const char* PDB = "pdb";
        const char* PRM = "prmtop,prm7";

        if (match(path, DMS)) return DmsFileFormat;
        if (match(path, MAE)) return MaeFileFormat;
        if (match(path, PDB)) return PdbFileFormat;
        if (match(path, PRM)) return ParmTopFileFormat;
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

    SystemPtr LoadWithFormat(std::string const& path, FileFormat format) {
        SystemPtr m;
        switch (format) {
            case DmsFileFormat: 
                m=ImportDMS(path); 
                break;
            case MaeFileFormat: 
                m=ImportMAE(path); 
                break;
            case PdbFileFormat: 
                m=ImportPDB(path); 
                break;
            case ParmTopFileFormat: 
                m=ImportPrmTop(path); 
                break;
            default:
                ;
        }
        return m;
    }

    SystemPtr Load(std::string const& path, FileFormat* opt_format) {
        FileFormat format = GuessFileFormat(path);
        if (opt_format) *opt_format = format;
        return LoadWithFormat(path, format);
    }

}}
