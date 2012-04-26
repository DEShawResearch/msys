#ifndef desres_msys_mae_vdwmap_hxx
#define desres_msys_mae_vdwmap_hxx

#include <fastjson/fastjson.hxx>
#include <map>
#include <string>
#include <vector>

namespace desres { namespace msys { namespace mae {

    /* vdw information comes from four different places in a ct block
     * (combining rule, vdwtypes, vdwtypes_combined, sites) and gets used 
     * to create at least two different nbody instances (nonbonded and pairs).
     * This class takes care of parsing the vdw parameters, and constructing 
     * a mapping between atoms and their vdw parameters.  
     */

    typedef std::vector<double> VdwParam;
    typedef std::string         VdwType;
    typedef std::pair<std::string,std::string>  CombinedKey;
    typedef std::map<CombinedKey, VdwParam>     CombinedMap;

    class VdwMap {

        std::string _funct;     /* ffio_funct in vdwtypes */
        std::string _rule;      /* combining rule */
        int         _nprops;    /* number of vdw properties */

        /* mapping from vdw typename to parameter set */
        typedef std::map<std::string, VdwParam> ParamMap;
        ParamMap    _params;

        /* vdw types for each site */
        typedef std::vector<std::string> VdwNameList;
        VdwNameList _vdwnames;
        VdwNameList _vdwnamesB;
        std::vector<double> _chargeB;

        /* combined params */
        CombinedMap _combined;

    public:
        explicit VdwMap( const fastjson::Json& ffio_ff );

        const std::string& funct() const { return _funct; }
        const std::string& rule() const  { return _rule;  }
        int nprops() const { return _nprops; }

        /* parameters for given name, containing the value from each
         * ffio_cN column found in the ffio_vdwtypes section */
        const VdwParam& param( const VdwType& name ) const;
        
        /* true if there is a combined param for the given pair of types */
        bool has_combined( const VdwType& name1, const VdwType& name2 ) const {
            return _combined.find(std::make_pair(name1,name2))!=_combined.end();
        }

        /* get the combined params for the given types */
        const VdwParam& param( const VdwType& name1, 
                               const VdwType& name2 ) const;

        /* return the combined params themselves so that they can be
         * exported. */
        const CombinedMap& combined() const { return _combined; }

        /* vdw type for given 0-based site id */
        const VdwType& type( int id ) const;

        /* alchemical vdw type, or empty string if not present */
        const VdwType& typeB( int id ) const;

        /* alchemical charge, or HUGE_VAL if not present */
        double chargeB( int id ) const;
    };

}}}

#endif
