#ifndef desres_msys_mae_ff_hxx
#define desres_msys_mae_ff_hxx

#include "sitemap.hxx"
#include "vdwmap.hxx"
#include "../term_table.hxx"
#include "../schema.hxx"
#include <stdexcept>
#include <sstream>
#include <boost/algorithm/string.hpp>

using namespace desres::msys;
using namespace desres::msys::mae;
using desres::fastjson::Json;

#define FFIO_ERROR( x ) do { \
    std::stringstream ss; \
    ss << x; \
    throw std::runtime_error(ss.str()); \
} while (0)

namespace desres { namespace msys { namespace mae { 

    using desres::fastjson::Json;

    struct Ffio {
        virtual ~Ffio() {}

        /* register an importer */
        static void put( const std::string& name, const Ffio * p );

        /* get importer for block, or NULL if none registered */
        static const Ffio* get( const std::string& name );

        /* if true, then the entire ffio_ff block is provided, not just
         * the block under which the importer is registered */
        virtual bool wants_all() const { return false; }

        virtual void apply( SystemPtr h,
                            const Json& blk, 
                            const SiteMap& sitemap,
                            const VdwMap&  vdwmap ) const = 0;
    };

    /* statically construct one of these for each plugin; e.g.,
     * namespace { RegisterFfio<VdwTypes> _("ffio_vdwtypes"); }
     */
    template <typename T>
    class RegisterFfio {
        Ffio * plugin;
    public:
        explicit RegisterFfio(const std::string& name ) {
            plugin = new T;
            Ffio::put( name, plugin );
        }
        ~RegisterFfio() { delete plugin; }
    };

    /* list of parameter values */
    typedef std::vector<double> ParamList;

    /* mae column names */
    typedef std::vector<std::string> StringList;

    /* mae columns */
    typedef std::vector<const Json *> ColumnList;

    /* This class supplies a hash on the actual parameter values for a
     * ParamTable, assuming they are all of FloatType. */
    class ParamMap {

        typedef std::map<ParamList, Id> ParamDict;

        ParamTablePtr   _params;        /* the tracked parameter table */
        ColumnList      _columns;       /* mae columns */
        ParamDict       _pdict;         /* map param values to param id */

    public:
        /* construct from param table, forcefield block, and names of the
         * columns to use for parameters.  If maecols is empty, then
         * "canonical" mae columns will be used, one for each property. */
        ParamMap( ParamTablePtr params, 
                  const Json& blk,
                  unsigned ncols,
                  const char** maecols );
        ParamMap( ParamTablePtr params, 
                  const Json& blk, 
                  bool alchemical=false);

        /* the parma table */
        ParamTablePtr params() const { return _params; }

        /* read the values from the given row in the forcefield block */
        Id add(int row);

        /* add a parameter with specified values not necessarily derived
         * from any row */
        Id add(const ParamList& param);
    };

    typedef boost::shared_ptr<ParamMap> ParamMapPtr;

}}}

#endif
