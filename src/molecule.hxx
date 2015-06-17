#ifndef desres_msys_molecule_hxx
#define desres_msys_molecule_hxx

#include <map>
#include <memory>
#include <string>

namespace desres { namespace msys {

    // a molecule is designed to efficiently represent small molecules, 
    // such as those found in sdfiles.  It does not support atom deletion,
    // bond deletion, etc.  
    class Molecule {
    public:
        struct Atom {
            char atomic_number = 0;
            char formal_charge = 0;
            char stereo_parity = 0;
            char aromatic      = 0;
            char implicit_h    = 0;
            float x = 0;
            float y = 0;
            float z = 0;
        };

        struct Bond {
            unsigned short i = 0;
            unsigned short j = 0;
            unsigned char order     = 0;
            unsigned char stereo    = 0;
            unsigned char aromatic  = 0;
        };

        typedef std::map<std::string,std::string> Data;

        unsigned short natoms() const { return _natoms; }
        unsigned short nbonds() const { return _nbonds; }

        const Atom& atom(int i) const { return _atoms[i]; }
        const Bond& bond(int i) const { return _bonds[i]; }
        Atom& atom(int i)             { return _atoms[i]; }
        Bond& bond(int i)             { return _bonds[i]; }

        typedef std::string Name;
        Name const& name() const { return _name; }
        Name&       name()       { return _name; }

        Data const& data() const { return _data; }
        Data&       data()       { return _data; }

        // construct new molecule - throw if size limits are exceeded
        Molecule(unsigned natoms, unsigned nbonds);

        // the usual boilerplate
        Molecule(Molecule const& c);
        Molecule& operator=(Molecule const& c);

    private:
        std::unique_ptr<Atom[]> _atoms;
        std::unique_ptr<Bond[]> _bonds;
        Data _data;
        Name _name;
        unsigned short _natoms;
        unsigned short _nbonds;
    };

    typedef std::unique_ptr<Molecule> MoleculePtr;
    class MoleculeIterator {
    public:
        virtual ~MoleculeIterator() {}
        virtual MoleculePtr next() = 0;
    };
    typedef std::unique_ptr<MoleculeIterator> MoleculeIteratorPtr;

}}

#endif
