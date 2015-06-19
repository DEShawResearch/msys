#ifndef desres_msys_smiles_smiles_hxx
#define desres_msys_smiles_smiles_hxx

#include "../molecule.hxx"
#include "../types.hxx"
#include <string>
#include <map>
#include <string.h>

namespace desres { namespace msys { namespace smiles {

    struct atom_t;
    struct chain_t;
    struct branch_t;
    struct ringbond_t;

    struct Smiles {
        const char* txt;
        int pos;
        inline char getc() { return txt[pos++]; }
        std::string error;

        /* mapping from rnum to <atom,bond> */
        typedef std::pair<atom_t*,char> rnum_t;
        typedef std::map<int,rnum_t> ringmap;
        ringmap rnums;

        void* scanner;
        std::vector<Molecule::Atom> atoms;
        std::vector<Molecule::Bond> bonds;

        Smiles();

        /* defined in smiles.l */
        void init_scanner();
        void destroy_scanner();

        MoleculePtr parse(std::string const& s);

        void addAtom(atom_t* a, bool organic);
        void add(atom_t* ai, atom_t* aj, char bond);
        void add(atom_t* a, branch_t* branches);
        void add(atom_t* a, ringbond_t* ringbonds);

        int addh(Id atm, int v1, int v2=0, int v3=0);
        void addh(Id atm);

        void finish(chain_t* chain);
    };

    struct ringbond_t {
        char bond;  /* 0 | '-' | '=' | '#' | '$' | ':' | '/' | '\'       */
        int  id;
        ringbond_t *next;

        ringbond_t(char b, int i) : bond(b), id(i), next() {}
    };

    struct branch_t {
        char bond;  /* 0 | '-' | '=' | '#' | '$' | ':' | '/' | '\' | '.' */
        chain_t *chain;
        branch_t *next;

        explicit branch_t(char b, chain_t* c) : bond(b), chain(c), next() {}
    };

    struct atom_t {
        atom_t *next;
        Id      id;                 /* System id once constructed */

        char          name[3];
        unsigned char isotope;
        unsigned char chiral;
        unsigned char hcount;
        signed char   charge;
        unsigned char klass;
        //char          bond;

        void setName(const char* s) { strcpy(name,s); }
        void clear() { memset(this,0,sizeof(*this)); }
        atom_t() { clear(); }
    };

    struct chain_t {
        atom_t *first;
        atom_t *last;

        explicit chain_t(atom_t* f) : first(f), last(f) {}
    };

}}}
#endif
