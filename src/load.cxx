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

}

namespace desres { namespace msys {

    SystemPtr Load(std::string const& path) {

        const char* DMS = "dms,dms.gz";
        const char* MAE = "mae,mae.gz,maegz,maeff,maeff.gz,cms,cms.gz";
        const char* PDB = "pdb";
        const char* PRM = "prmtop,prm7";

        if (match(path, DMS)) return ImportDMS(path);
        if (match(path, MAE)) return ImportMAE(path);
        if (match(path, PDB)) return ImportPDB(path);
        if (match(path, PRM)) return ImportPrmTop(path);
        MSYS_FAIL("Could not guess file type of '" << path << "'");
    }

}}
