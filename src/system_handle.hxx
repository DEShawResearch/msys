#ifndef desres_msys_system_handle_hxx
#define desres_msys_system_handle_hxx

#include "atom.hxx"
#include "bond.hxx"
#include "residue.hxx"
#include "chain.hxx"
#include "param_handle.hxx"
#include "term_handle.hxx"

namespace desres { namespace msys {

    class SystemHandle {
    
        SystemPtr _system;
    
    public:
        /* default constructor creates invalid handle */
        SystemHandle() {}
        SystemHandle(SystemPtr system) : _system(system) {}

        SystemPtr ptr() const { return _system; }

        bool operator==(const SystemHandle& sys) const {
            return _system==sys._system;
        }

        GlobalCell& cell() {
            return _system->global_cell;
        }

        NonbondedInfo nonbondedInfo() const {
            return _system->nonbonded_info;
        }
        void setNonbondedInfo(const NonbondedInfo& info) {
            _system->nonbonded_info = info;
        }

        /* accessors for full set of system elements */
        AtomList atoms();
        BondList bonds();
        ResidueList residues();
        ChainList chains();

        /* delete a list of elements */
        void delAtoms(const AtomList& atoms);
        void delBonds(const BondList& bonds);
        void delResidues(const ResidueList& residues);
        void delChains(const ChainList& chains);

        /* accessors for individual structural elements */
        Atom atom(Id id);
        Bond bond(Id id);
        Residue residue(Id id);
        Chain chain(Id id);

        /* add elements to the system, creating parent entities as needed */
        Atom addAtom();
        Residue addResidue();
        Chain addChain(); 

        /* counts */
        Id atomCount() { return _system->atomCount(); }
        Id bondCount() { return _system->bondCount(); }
        Id residueCount() { return _system->residueCount(); }
        Id chainCount() { return _system->chainCount(); }

        /* atom properties */
        Id addAtomProp(const String& name, ValueType type) {
            return _system->atomProps()->addProp(name, type);
        }
        Id atomPropIndex(const String& name) const {
            return _system->atomProps()->propIndex(name);
        }
        bool hasAtomProp(const String& name) const {
            return atomPropIndex(name)!=BadId;
        }
        Id atomPropCount() const {
            return _system->atomProps()->propCount();
        }
        String atomPropName(Id i) const {
            return _system->atomProps()->propName(i);
        }
        ValueType atomPropType(Id i) const {
            return _system->atomProps()->propType(i);
        }
        ValueRef atomProp(Id id, Id col) const {
            return _system->atomProps()->value(id,col);
        }
        /* The const-ness of this method is a big lie. */
        ValueRef atomProp(Id id, const String& name) const {
            return atomProp(id, atomPropIndex(name));
        }

        /* add a bond between existing atoms */
        Bond addBond(Atom a1, Atom a2);

        /* find the bond between the given atoms; return invalid Bond if none */
        Bond findBond(Atom a1, Atom a2);

        /* bond properties */
        Id addBondProp(const String& name, ValueType type) {
            return _system->bondProps()->addProp(name, type);
        }
        Id bondPropIndex(const String& name) const {
            return _system->bondProps()->propIndex(name);
        }
        bool hasBondProp(const String& name) const {
            return bondPropIndex(name)!=BadId;
        }
        Id bondPropCount() const {
            return _system->bondProps()->propCount();
        }
        String bondPropName(Id i) const {
            return _system->bondProps()->propName(i);
        }
        ValueType bondPropType(Id i) const {
            return _system->bondProps()->propType(i);
        }
        ValueRef bondProp(Id id, Id col) const {
            return _system->bondProps()->value(id,col);
        }
        ValueRef bondProp(Id id, const String& name) const {
            return bondProp(id, bondPropIndex(name));
        }

        TermTableHandle addTable(const String& name, Id natoms,
                              ParamTableHandle params = ParamTableHandle()) {
            return _system->addTable(name,natoms,params.ptr());
        }
        std::vector<String> tableNames() const {
            return _system->tableNames();
        }
        TermTableHandle table(const String& name) const {
            return _system->table(name);
        }
        bool hasTable(const String& name) const {
            return _system->table(name);
        }
        void delTable(const String& name) {
            _system->delTable(name);
        }
        void removeTable(TermTableHandle table) {
            _system->removeTable(table.ptr());
        }

        AtomList atomselect( String const& seltext ) const;

        AtomList appendSystem( SystemHandle src );

        TermTableHandle addTableFromSchema( String const& type,
                                            String const& name="");

        TermTableHandle addNonbondedFromSchema( String const& funct,
                                                String const& rule );

    };

    SystemHandle CreateSystem();
    SystemHandle CloneSystem( AtomList const& atoms );

    SystemHandle LoadMAE(String const& path,
                           bool ignore_unrecognized = false );
    SystemHandle LoadDMS(String const& path, 
                           bool with_forcefield = true );

    void SaveDMS(SystemHandle h, String const& path);
    void SaveMAE(SystemHandle h, String const& path);


    template <class Obj>
    SystemHandle system(Obj& obj) {
        return SystemHandle(obj.sys()->shared_from_this());
    }

}}

#endif
