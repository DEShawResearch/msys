#ifndef desres_msys_system_hxx
#define desres_msys_system_hxx

#include <vector>
#include <map>
#include <boost/enable_shared_from_this.hpp>

#include "types.hxx"
#include "term_table.hxx"
#include "param_table.hxx"

namespace desres { namespace msys {

    class TermTable;
    typedef boost::shared_ptr<TermTable> TermTablePtr;

        
    struct Vec3 {
        Float x,y,z;
        
        Vec3() 
        : x(0), y(0), z(0) {}

        Vec3(Float const& _x, Float const& _y, Float const& _z)
        : x(_x), y(_y), z(_z) {}

        Float& operator[](unsigned i) {
            return i==0 ? x : i==1 ? y : z;
        }
        Float const& operator[](unsigned i) const {
            return i==0 ? x : i==1 ? y : z;
        }

        bool operator==(Vec3 const& o) const {
            return x==o.x && y==o.y && z==o.z;
        }
        bool operator!=(Vec3 const& o) const {
            return x!=o.x || y!=o.y || z!=o.z;
        }

    };

    struct GlobalCell {
        Vec3 A, B, C;
        Vec3& operator[](unsigned i) { 
          return i==0 ? A : i==1 ? B : C;
        }
        Vec3 const& operator[](unsigned i) const { 
          return i==0 ? A : i==1 ? B : C;
        }
    };

    struct NonbondedInfo {
        String vdw_funct;
        String vdw_rule;
    };

    struct atom_t {
        Id  gid;        /* For dms, the primary key of the particle */
        Id  fragid;
        Id  residue;
    
        Float x,y,z;    /* position */
        Float charge;   /* partial charge */
        Float vx,vy,vz; /* velocity */
        Float mass;
        Float chargeB;
    
        int atomic_number;
        int formal_charge;
        int moiety;
        bool alchemical;
    
        String name;
    
        atom_t() 
        : gid(BadId), fragid(BadId), residue(BadId),
          x(0), y(0), z(0), charge(0), vx(0), vy(0), vz(0), mass(0), chargeB(0),
          atomic_number(0), formal_charge(0), moiety(0), alchemical(false)
        {}
    };
    
    struct bond_t {
        Id  i;      /* id of first bond partner     */
        Id  j;      /* id of second bond partner    */
        Float order;  /* resonant bond order          */
    
        bond_t() : i(BadId), j(BadId), order(1) {}
        bond_t(Id ai, Id aj) : i(ai), j(aj), order(1) {}
        Id other(Id id) const { return id==i ? j : i; }
    };
    
    struct residue_t {
        String  name;
        int     num;
        Id      chain;
    
        residue_t() : num(0), chain(BadId) {}
    };
    
    struct chain_t {
        String  name;
    
        chain_t() {}
    };
    
    class System : public boost::enable_shared_from_this<System> {
    
        static IdList _empty;
    
        /* _atoms maps an id to an atom.  We almost never delete atoms, so
         * keep track of deleted atoms in a separate set.  This is needed only
         * by an atom iterator */
        typedef std::vector<atom_t> AtomList;
        AtomList    _atoms;
        IdSet       _deadatoms;

        /* additional properties for atoms */
        ParamTablePtr _atomprops;

        /* same deal for bonds */
        typedef std::vector<bond_t> BondList;
        BondList    _bonds;
        IdSet       _deadbonds;
        ParamTablePtr _bondprops;
    
        /* map from atom id to 0 or more bond ids.  We do keep this updated when
         * atoms or bonds are deleted */
        MultiIdList   _bondindex;
    
        typedef std::vector<residue_t> ResidueList;
        ResidueList _residues;
        IdSet       _deadresidues;
        MultiIdList   _residueatoms;  /* residue id -> atom ids */
    
        typedef std::vector<chain_t> ChainList;
        ChainList   _chains;
        IdSet       _deadchains;
        MultiIdList   _chainresidues; /* chain id -> residue id */
    
        typedef std::map<String,TermTablePtr> TableMap;
        TableMap    _tables;

        /* extra tables.  basically a hack for cmap */
        typedef std::map<String, ParamTablePtr> ExtraMap;
        ExtraMap    _extras;

        /* create only as shared pointer. */
        System();
    public:
        static boost::shared_ptr<System> create();
        ~System();

        String          name;
        GlobalCell      global_cell;
        NonbondedInfo   nonbonded_info;
    
        /* element accessors */
        atom_t& atom(Id id) { return _atoms.at(id); }
        bond_t& bond(Id id) { return _bonds.at(id); }
        residue_t& residue(Id id) { return _residues.at(id); }
        chain_t& chain(Id id) { return _chains.at(id); }
    
        const atom_t& atom(Id id) const { return _atoms.at(id); }
        const bond_t& bond(Id id) const { return _bonds.at(id); }
        const residue_t& residue(Id id) const { return _residues.at(id); }
        const chain_t& chain(Id id) const { return _chains.at(id); }

        /* add an element */
        Id addAtom(Id residue);
        Id addBond(Id i, Id j);
        Id addResidue(Id chain);
        Id addChain();

        /* delete an element */
        void delAtom(Id id);
        void delBond(Id id);
        void delResidue(Id id);
        void delChain(Id id);

        /* One more than highest valid id */
        Id maxAtomId() const;
        Id maxBondId() const;
        Id maxResidueId() const;
        Id maxChainId() const;

        /* list of elements ids */
        IdList atoms() const;
        IdList bonds() const;
        IdList residues() const;
        IdList chains() const;

        /* number of elements */
        Id atomCount() const { return _atoms.size() - _deadatoms.size(); }
        Id bondCount() const { return _bonds.size() - _deadbonds.size(); }
        Id residueCount() const {return _residues.size()-_deadresidues.size();}
        Id chainCount() const { return _chains.size() - _deadchains.size(); }

        /* count of subelements */
        Id bondCountForAtom(Id id) const {
            if (id>=_bondindex.size()) return 0;
            return _bondindex[id].size();
        }
        Id atomCountForResidue(Id id) const {
            if (id>=_residueatoms.size()) return 0;
            return _residueatoms[id].size();
        }
        Id residueCountForChain(Id id) const {
            if (id>=_chainresidues.size()) return 0;
            return _chainresidues[id].size();
        }

        /* list of subelements */
        IdList const& bondsForAtom(Id id) const {
            if (id>=_bondindex.size()) return _empty;
            return _bondindex[id];
        }
        IdList const& atomsForResidue(Id id) const {
            if (id>=_residueatoms.size()) return _empty;
            return _residueatoms[id];
        }
        IdList const& residuesForChain(Id id) const {
            if (id>=_chainresidues.size()) return _empty;
            return _chainresidues[id];
        }

        /* is the the given element id valid? */
        bool hasAtom(Id id) const {
            return id<_atoms.size() && !_deadatoms.count(id);
        }
        bool hasBond( Id id ) const {
            return id<_bonds.size() && !_deadbonds.count(id);
        }
        bool hasResidue(Id id) const {
            return id<_residues.size() && !_deadresidues.count(id);
        }
        bool hasChain(Id id) const {
            return id<_chains.size() && !_deadchains.count(id);
        }


        /* delete multiple elements at once (convenience) */
        template <typename T> void delAtoms(const T& ids) {
            typedef typename T::const_iterator iterator;
            for (iterator i=ids.begin(), e=ids.end(); i!=e; ++i) delAtom(*i);
        }
        template <typename T> void delBonds(const T& ids) {
            typedef typename T::const_iterator iterator;
            for (iterator i=ids.begin(), e=ids.end(); i!=e; ++i) delBond(*i);
        }
        template <typename T> void delResidues(const T& ids) {
            typedef typename T::const_iterator iterator;
            for (iterator i=ids.begin(), e=ids.end(); i!=e; ++i) delResidue(*i);
        }
        template <typename T> void delChains(const T& ids) {
            typedef typename T::const_iterator iterator;
            for (iterator i=ids.begin(), e=ids.end(); i!=e; ++i) delChain(*i);
        }

        /* operations on term tables */
        std::vector<String> tableNames() const;
        /* fetch the table with the given name; return NULL if not present */
        TermTablePtr table(const String& name) const;
        /* get the name of the table; throw if table doesn't belong to this */
        String tableName(boost::shared_ptr<TermTable const> table) const;
        TermTablePtr addTable(const String& name, Id natoms,
                              ParamTablePtr ptr = ParamTablePtr() );
        void delTable(const String& name);
        void removeTable(TermTablePtr terms);

        /* operations on extra tables */
        std::vector<String> extraNames() const;
        ParamTablePtr extra(const String& name) const;
        void addExtra(const String& name, ParamTablePtr ptr);
        void delExtra(const String& name);
        void removeExtra(ParamTablePtr extra);

        /* assign the atom to the given residue */
        void setResidue(Id atom, Id residue);

        /* extended atom properties */
        Id atomPropCount() const;
        String atomPropName(Id i) const;
        ValueType atomPropType(Id i) const;
        Id atomPropIndex(String const& name) const;
        Id addAtomProp(String const& name, ValueType type);
        void delAtomProp(Id index);
        ValueRef atomPropValue(Id term, Id index);
        ValueRef atomPropValue(Id term, String const& name);

        /* extended bond properties */
        Id bondPropCount() const;
        String bondPropName(Id i) const;
        ValueType bondPropType(Id i) const;
        Id bondPropIndex(String const& name) const;
        Id addBondProp(String const& name, ValueType type);
        void delBondProp(Id index);
        ValueRef bondPropValue(Id term, Id index);
        ValueRef bondPropValue(Id term, String const& name);

        /* find a bond between i and j (order independent), returning
         * its id or BadId if such a bond does not exist. */
        Id findBond( Id i, Id j) const;

        /* ids of atoms bonded to given atom */
        IdList bondedAtoms(Id id) const;
    
        /* update the fragid of each atom according to its bond topology:
         * bonded atoms share the same fragid.  Return the number of
         * frags found, and atomid to fragment partitioning if requested */
        Id updateFragids(MultiIdList* fragments=NULL);
    };

    typedef boost::shared_ptr<System> SystemPtr;

}}

#endif
