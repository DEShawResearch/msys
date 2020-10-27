#ifndef desres_msys_annotated_system_hxx
#define desres_msys_annotated_system_hxx

#include "system.hxx"

namespace desres { namespace msys {

    class AnnotatedSystem {
    public:
        enum Flags { Default            = 0 
                   , AllowBadCharges    = 1 << 0
        };

    private:
        struct atom_data_t {
            static const unsigned char MaxBonds = 6;

            bool          aromatic = false;
            unsigned char hcount = 0;
            unsigned char valence = 0;
            unsigned char degree = 0;
            unsigned char lone_electrons = 0;
            unsigned char ring_bonds = 0;
            unsigned char hybridization = 0;
            unsigned char atomic_number = 0;
            char formal_charge = 0;

            Id bond[MaxBonds];
            IdList rings_idx;
        };

        struct bond_data_t {
            Id i=BadId;
            Id j=BadId;
            bool aromatic = false;
            unsigned char order = 0;
            IdList rings_idx;

            Id other(Id id) const { return id==i ? j : i; }
        };
        struct ring_t {
            IdList atoms;
            IdList bonds;
        };
        struct ring_system_t {
            IdList atoms;
            IdList bonds;
            IdList rings;
        };

        std::vector<atom_data_t> _atoms;
        std::vector<bond_data_t> _bonds;
        std::vector<ring_t> _rings;
        std::vector<ring_system_t> _ring_systems;

        std::vector<String> _errors;

        /* Helper functions for constructor */
        void compute_ring_systems(SystemPtr sys);
        bool is_aromatic(SystemPtr, const IdList& atms, const IdList& bnds);
        void compute_aromaticity(SystemPtr);

    public:
        /* Create an annotated system. sys must have correct bond orders
         * and formal charges (assigned using AssignBondOrderAndFormalCharge
         * or otherwise) */
        AnnotatedSystem(SystemPtr sys, unsigned flags = Default);

        /* make non-copyable */
        AnnotatedSystem(AnnotatedSystem const&) = delete;
        AnnotatedSystem& operator=(AnnotatedSystem const&) = delete;
        ~AnnotatedSystem() = default;

        std::vector<std::string> errors() const { return _errors; }

        // non-deleted, non-pseudo atom ids
        IdList atoms() const;
        Id atomCount() const { return _atoms.size(); }
        Id bondCount() const { return _bonds.size(); }
        atom_data_t const& atomFAST(Id atom) const {
            return _atoms[atom];
        }
        bond_data_t const& bondFAST(Id bond) const {
            return _bonds[bond];
        }

        /* Is atom aromatic */
        bool atomAromatic(Id atom) const {
            return _atoms.at(atom).aromatic; }
        /* Number of attached hydrogens */
        int atomHcount(Id atom) const {
            return _atoms.at(atom).hcount; }
        /* Sum of bond orders of non-pseudo bonds */
        int atomValence(Id atom) const {
            return _atoms.at(atom).valence; }
        /* Total number of non-pseudo bonds */
        int atomDegree(Id atom) const {
            return _atoms.at(atom).degree; }
        int atomLoneElectrons(Id atom) const {
            return _atoms.at(atom).lone_electrons; }
        /* Hybridization (1 = sp, 2 = sp2, 3 = sp3, 4 = sp3d, etc.) */
        int atomHybridization(Id atom) const {
            return _atoms.at(atom).hybridization; }
        /* Total number of ring bonds */
        int atomRingBonds(Id atom) const {
            return _atoms.at(atom).ring_bonds; }
        /* Total number of rings */
        int atomRingCount(Id atom) const {
            return _atoms.at(atom).rings_idx.size(); }
        /* Is atom in ring of particular size */
        bool atomInRingSize(Id atom, unsigned size) const {
            for (Id ring : _atoms.at(atom).rings_idx)
                if (_rings.at(ring).atoms.size() == size) return true;
            return false;
        }
        bool atomInSmallestRingSize(Id atom, unsigned size) const {
            size_t ring_size = BadId;
            for (Id ring : _atoms.at(atom).rings_idx)
                ring_size = std::min(ring_size, _rings.at(ring).atoms.size());
            return ring_size == size;
        }
        /* SSSR rings containing this atom */
        MultiIdList atomRings(Id atom) const {
            MultiIdList rings;
            for (Id ring : _atoms.at(atom).rings_idx)
                rings.push_back(_rings.at(ring).atoms);
            return rings;
        }

        /* Is bond aromatic */
        bool bondAromatic(Id bond) const {
            return _bonds.at(bond).aromatic; }
        /* Total number of rings */
        int bondRingCount(Id bond) const {
            return _bonds.at(bond).rings_idx.size(); }
        bool bondInRingSize(Id bond, unsigned size) const {
            for (Id ring : _bonds.at(bond).rings_idx)
                if (_rings.at(ring).atoms.size() == size) return true;
            return false;
        }
        /* SSSR rings containing this bond */
        MultiIdList bondRings(Id bond) const {
            MultiIdList rings;
            for (Id ring : _bonds.at(bond).rings_idx)
                rings.push_back(_rings.at(ring).atoms);
            return rings;
        }

        /* All SSSR rings */
        int ringCount() const {
            return _rings.size();
        }
        MultiIdList rings() const {
            MultiIdList rings;
            for (const ring_t& ring : _rings)
                rings.push_back(ring.atoms);
            return rings;
        }
    };
}}

#endif
