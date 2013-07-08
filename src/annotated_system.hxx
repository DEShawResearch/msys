#ifndef desres_msys_annotated_system_hxx
#define desres_msys_annotated_system_hxx

#include "system.hxx"

namespace desres { namespace msys {

    class AnnotatedSystem :
        public boost::enable_shared_from_this<AnnotatedSystem> {

        private:
            struct atom_data_t {
                bool aromatic;
                unsigned char hcount;
                unsigned char valence;
                unsigned char degree;
                unsigned char lone_pairs;
                unsigned char ring_bonds;
                unsigned char hybridization;
                IdList rings_idx;

                atom_data_t() 
                : aromatic(), hcount(), valence(), degree(),
                  lone_pairs(), ring_bonds(), hybridization()
                {}
            };
            struct bond_data_t {
                bool aromatic;
                IdList rings_idx;

                bond_data_t()
                : aromatic()
                {}
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

            SystemPtr _sys;
            std::vector<atom_data_t> _atoms;
            std::vector<bond_data_t> _bonds;
            std::vector<ring_t> _rings;
            std::vector<ring_system_t> _ring_systems;

            AnnotatedSystem(SystemPtr sys);

            /* Helper functions for constructor */
            void compute_ring_systems();
            bool is_aromatic(const IdList& atoms, const IdList& bonds);
            void compute_aromaticity();

        public:
            /* Create an annotated system. sys must have correct bond orders
             * and formal charges (assigned using AssignBondOrderAndFormalCharge
             * or otherwise) */
            static boost::shared_ptr<AnnotatedSystem> create(SystemPtr sys) {
                return boost::shared_ptr<AnnotatedSystem>(
                        new AnnotatedSystem(sys));
            }
            SystemPtr system() const { return _sys; }

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
            int atomLonePairs(Id atom) const {
                return _atoms.at(atom).lone_pairs; }
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
                BOOST_FOREACH(Id ring, _atoms.at(atom).rings_idx)
                    if (_rings.at(ring).atoms.size() == size) return true;
                return false;
            }
            /* SSSR rings containing this atom */
            void atomRings(Id atom, MultiIdList& rings) const {
                rings.clear();
                BOOST_FOREACH(Id ring, _atoms.at(atom).rings_idx)
                    rings.push_back(_rings.at(ring).atoms);
            }

            /* Is bond aromatic */
            bool bondAromatic(Id bond) const {
                return _bonds.at(bond).aromatic; }
            /* Total number of rings */
            int bondRingCount(Id bond) const {
                return _bonds.at(bond).rings_idx.size(); }
            bool bondInRingSize(Id bond, unsigned size) const {
                BOOST_FOREACH(Id ring, _bonds.at(bond).rings_idx)
                    if (_rings.at(ring).atoms.size() == size) return true;
                return false;
            }
            /* SSSR rings containing this bond */
            void bondRings(Id bond, MultiIdList& rings) const {
                rings.clear();
                BOOST_FOREACH(Id ring, _bonds.at(bond).rings_idx)
                    rings.push_back(_rings.at(ring).atoms);
            }

            /* All SSSR rings */
            int ringCount() const {
                return _rings.size(); }
            void rings(MultiIdList& rings) const {
                rings.clear();
                BOOST_FOREACH(const ring_t& ring, _rings)
                    rings.push_back(ring.atoms);
            }
    };
    typedef boost::shared_ptr<AnnotatedSystem> AnnotatedSystemPtr;
}}

#endif
