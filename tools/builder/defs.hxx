#ifndef desres_msys_builder_defs_hxx
#define desres_msys_builder_defs_hxx

#include <string>
#include <vector>
#include <map>
#include <iostream>
#include "../../src/types.hxx"

namespace desres { namespace msys { namespace builder {

    struct type_t {
        int     id;
        int     anum;
        double  mass;

        type_t() : id(), anum(), mass() {}
    };
    typedef std::map<std::string, type_t> TypeMap;

    struct adef_t {
        std::string name;
        int         res;
        int         rel;

        adef_t() : res(), rel() {}

        /* if check_res is true, then interpret a leading digit as a
         * residue specification for a patch residue; otherwise include
         * it in the atom name. */
        void parse(const char *s, bool check_res = true);
    };

    struct atom_t {
        adef_t      def;
        std::string type;
        double      charge;

        atom_t() : charge() {}
    };

    struct bond_t {
        adef_t      def1;
        adef_t      def2;

        bond_t() {}
    };

    struct conf_t {
        adef_t  def1;
        adef_t  def2;
        adef_t  def3;
        adef_t  def4;

        double  r1;
        double  a1;
        double  phi;
        double  a2;
        double  r2;

        bool    improper;

        conf_t() : r1(), a1(), phi(), a2(), r2(), improper() {}
    };

    struct resdef_t {
        std::vector<atom_t> atoms;
        std::vector<bond_t> bonds;
        std::vector<conf_t> confs;

        std::vector<adef_t> delatoms;
        std::vector<bond_t> delbonds;

        std::string         pfirst;
        std::string         plast;
        
        bool                patch;

        resdef_t() : patch() {}

        void patch_topology(resdef_t& topo) const;

        Id atom_index(std::string const& name) const {
            for (Id i=0; i<atoms.size(); i++) {
                if (atoms[i].def.name==name) return i;
            }
            return BadId;
        }
    };

    typedef std::map<std::string, resdef_t> ResdefMap;

    struct file_t {
        std::string     path;
    };
    typedef std::vector<file_t> FileList;

    struct defs_t {
        std::string     pfirst;
        std::string     plast;

        FileList        files;
        TypeMap         types;
        ResdefMap       resdefs;

        defs_t() {}

        resdef_t& add_resdef(std::string const& name) {
            return resdefs[name];
        }

        type_t& add_type(std::string const& name) {
            return types[name];
        }

        void import_charmm_topology(std::string const& path);
    };

}}}
#endif
