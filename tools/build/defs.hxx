#ifndef desres_msys_builder_defs_hxx
#define desres_msys_builder_defs_hxx

#include <string>
#include <vector>
#include <map>
#include <iostream>

namespace desres { namespace msys { namespace builder {

    struct type_t {
        std::string element;
        int         id;
        double      mass;

        type_t() : id(), mass() {}
    };
    typedef std::map<std::string, type_t> TypeMap;

    struct adef_t {
        std::string name;
        int         res;
        int         rel;

        adef_t() : res(), rel() {}

        void parse(const char *s);
    };

    struct atom_t {
        adef_t      def;
        bool        del;
        std::string type;
        double      charge;

        atom_t() : charge() {}
    };

    struct bond_t {
        adef_t      def1;
        adef_t      def2;
        bool        del;

        bond_t() : del() {}
    };

    struct conf_t {
        adef_t  def1;
        adef_t  def2;
        adef_t  def3;
        adef_t  def4;
        bool    del;

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

        std::string         pfirst;
        std::string         plast;
        
        bool        patch;

        resdef_t() : patch() {}

        atom_t& add_atom() {
            atoms.resize(atoms.size()+1);
            return atoms.back();
        }
        bond_t& add_bond() {
            bonds.resize(bonds.size()+1);
            return bonds.back();
        }
        conf_t& add_conf() {
            confs.resize(confs.size()+1);
            return confs.back();
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
