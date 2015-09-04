#ifndef desres_msys_system_hxx
#define desres_msys_system_hxx

#include <vector>
#include <map>

#include "term_table.hxx"
#include "provenance.hxx"
#include "value.hxx"
#include "smallstring.hxx"
#include "pfx/graph.hxx"

#include <boost/serialization/binary_object.hpp>

namespace desres { namespace msys {

    class TermTable;
    typedef boost::shared_ptr<TermTable> TermTablePtr;

    class ParamTable;
    typedef boost::shared_ptr<ParamTable> ParamTablePtr;

    class GlobalCell {
        std::vector<double> data;

        friend class boost::serialization::access;
        template <typename Ar> void serialize(Ar& a, unsigned) { a & data; }
    public:
        GlobalCell() : data(9) {}
        double*       operator[](unsigned i)       { return &data.at(3*i); }
        const double* operator[](unsigned i) const { return &data.at(3*i); }


        /* Allow merge when corresponding vectors are equal, or when
         * self or other is identically 0, in which case the non-zero
         * cell is adopted. */
        void merge(GlobalCell const& other);
    };

    struct NonbondedInfo {
        String vdw_funct;
        String vdw_rule;
        String es_funct;

        /* Allow merge when corresponding fields are equal, or when one
         * of them is empty, in which case the non-empty value is adopted.
         */
        void merge(NonbondedInfo const& other);

        template <typename Ar> void serialize(Ar& a, unsigned) {
            a & vdw_funct & vdw_rule & es_funct;
        }
    };

    enum AtomType {
        AtomOther   = 0,
        AtomProBack = 1,
        AtomNucBack = 2,
        AtomProSide = 3
    };

    struct atom_t {
        Id  fragid;
        Id  residue;
        int atomic_number;
        int formal_charge;
    
        Float x,y,z;    /* position */
        Float charge;   /* partial charge */
        Float vx,vy,vz; /* velocity */
        Float mass;
    
        SmallString<30> name;
        AtomType type;
        Float resonant_charge;
    
        atom_t() 
        : fragid(BadId), residue(BadId), atomic_number(0), formal_charge(0),
          x(0), y(0), z(0), charge(0), vx(0), vy(0), vz(0), mass(0), type(),
          resonant_charge()
        {}

        template <class Ar>
        void serialize(Ar& a, const unsigned int) {
            a & fragid & residue & atomic_number & formal_charge
              & x & y & z & charge 
              & vx & vy & vz & mass
              & name & type & resonant_charge
              ;
        }

        /* Don't abuse these.  In particular, bear in mind that that an atom's
         * memory location will move around if atoms are added.  */
        double       *pos()       { return &x; }
        const double *pos() const { return &x; }
    };
    
    struct bond_t {
        Id  i;                  /* id of first bond partner     */
        Id  j;                  /* id of second bond partner    */
        int order;              /* formal bond order            */
        Float resonant_order;   /* resonant bond order          */
    
        bond_t() : i(BadId), j(BadId), order(1), resonant_order(1) {}
        bond_t(Id ai, Id aj) : i(ai), j(aj), order(1), resonant_order(1) {}
        Id other(Id id) const { return id==i ? j : i; }

        template <class Ar>
        void serialize(Ar& a, const unsigned int) {
            a & i & j & order & resonant_order;
        }
    };
    
    enum ResidueType {
        ResidueOther    = 0,
        ResidueProtein  = 1,
        ResidueNucleic  = 2,
        ResidueWater    = 3
    };

    struct residue_t {
        Id      chain;
        int     resid;
        SmallString<30> name;
        SmallString<6> insertion;
        ResidueType type;
    
        residue_t() : chain(BadId), resid(), type() {}
    };
    
    struct chain_t {
        Id      ct;
        String  name;
        String  segid;
    
        chain_t() : ct(BadId) {}

        template <class Ar>
        void serialize(Ar& a, const unsigned int) {
            a & ct & name & segid;
        }
    };

    class component_t {
        ParamTablePtr _kv;

    public:
        /* constructor: maintain a single row with msys_name as the first
         * property. */
        component_t();

        /* must implement copy constructor so that we don't share _kv! */
        component_t(component_t const& c) { *this = c; }
        component_t& operator=(component_t const& c);

        /* getter/setter for name */
        String name() const;
        void setName(String const& s);

        /* other keys besides name */
        std::vector<String> keys() const;
        void del(String const& key);
        Id add(String const& key, ValueType type);
        ValueType type(String const& key) const;

        bool has(String const& key) const;
        ValueRef value(String const& key);
        ValueRef value(Id key);

        template <class Ar>
        void serialize(Ar& a, const unsigned int) {
            a & _kv;
        }
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
    
        typedef std::vector<component_t> CtList;
        CtList      _cts;
        IdSet       _deadcts;
        MultiIdList _ctchains; /* ct id -> chain id */
    
        typedef std::map<String,TermTablePtr> TableMap;
        TableMap    _tables;

        /* auxiliary tables.  basically a hack for cmap */
        typedef std::map<String, ParamTablePtr> AuxTableMap;
        AuxTableMap _auxtables;

        /* provenance.  Ideally, you would append to this just before 
         * serializing to disk. */
        std::vector<Provenance> _provenance;

        /* create only as shared pointer. */
        System();

        /* Cerealization */
        friend class boost::serialization::access;
        template <typename Ar> void serialize(Ar& a, unsigned) {
            using boost::serialization::binary_object;
            a & _deadatoms & _atomprops;
            Id an, bn, rn;
            if (Ar::is_saving::value) {
                an = _atoms.size();
                bn = _bonds.size();
                rn = _residues.size();
                a & an & bn & rn;
            } else {
                a & an & bn & rn;
                _atoms.resize(an);
                _bonds.resize(bn);
                _residues.resize(rn);
            }
            if (an) {
                binary_object b(&_atoms[0],an*sizeof(atom_t));
                a & b;
            }
            if (bn) {
                binary_object b(&_bonds[0],bn*sizeof(bond_t));
                a & b;
            }
            if (rn) {
                binary_object b(&_residues[0],rn*sizeof(residue_t));
                a & b;
            }

            a & _deadbonds & _bondprops;
            a & _deadresidues;
            a & _chains & _deadchains;
            a & _cts & _deadcts;
            a & _tables & _auxtables;
            a & _provenance;
            a & name & global_cell & nonbonded_info;
            if (Ar::is_loading::value) restore_indices();
        }
        void restore_indices();

    public:
        static boost::shared_ptr<System> create();
        ~System();

        String          name;
        GlobalCell      global_cell;
        NonbondedInfo   nonbonded_info;

        /* get the provenance history */
        std::vector<Provenance> const& provenance() const {
            return _provenance;
        }

        /* add a provenance entry.  Should be used only by ImportDMS! */
        void addProvenance(Provenance const& p) {
            _provenance.push_back(p);
        }

        /* element accessors */
        atom_t& atom(Id id) { return _atoms.at(id); }
        bond_t& bond(Id id) { return _bonds.at(id); }
        residue_t& residue(Id id) { return _residues.at(id); }
        chain_t& chain(Id id) { return _chains.at(id); }
        component_t& ct(Id id) { return _cts.at(id); }
    
        const atom_t& atom(Id id) const { return _atoms.at(id); }
        const bond_t& bond(Id id) const { return _bonds.at(id); }
        const residue_t& residue(Id id) const { return _residues.at(id); }
        const chain_t& chain(Id id) const { return _chains.at(id); }
        const component_t& ct(Id id) const { return _cts.at(id); }

        /* unchecked element accessors */
        atom_t& atomFAST(Id id) { return _atoms[id]; }
        bond_t& bondFAST(Id id) { return _bonds[id]; }
        residue_t& residueFAST(Id id) { return _residues[id]; }
        chain_t& chainFAST(Id id) { return _chains[id]; }
        component_t& ctFAST(Id id) { return _cts[id]; }

        const atom_t& atomFAST(Id id) const { return _atoms[id]; }
        const bond_t& bondFAST(Id id) const { return _bonds[id]; }
        const residue_t& residueFAST(Id id) const { return _residues[id]; }
        const chain_t& chainFAST(Id id) const { return _chains[id]; }
        const component_t& ctFAST(Id id) const { return _cts[id]; }

        /* id iterator, skipping deleted ids */
        class iterator {
            friend class System;
            Id            _i;
            const IdSet*  _dead;

            iterator(Id i, const IdSet* dead) 
            : _i(i), _dead(dead) {
                while (_dead && _dead->count(_i)) ++_i;
            }


            bool equal(iterator const& c) const { return _i==c._i; }
            const Id& dereference() const { return _i; }
            void increment() { do ++_i; while (_dead && _dead->count(_i)); }

        public:
            typedef std::forward_iterator_tag iterator_category;
            typedef Id value_type;
            typedef ptrdiff_t difference_type;
            typedef const Id* pointer;
            typedef const Id& reference;

            iterator() : _i(), _dead() {}
            const Id& operator*() const { return dereference(); }
            const Id* operator->() const { return &dereference(); }
            iterator& operator++() { increment(); return *this; }
            bool operator==(iterator const& c) const { return equal(c); }
            bool operator!=(iterator const& c) const { return !equal(c); }
        };

        template <typename T>
        void getPositions(T setter) const {
            for (iterator i=atomBegin(), e=atomEnd(); i!=e; ++i) {
                atom_t const& a = atomFAST(*i);
                *setter++ = a.x;
                *setter++ = a.y;
                *setter++ = a.z;
            }
        }

        template <typename T>
        void setPositions(T getter) {
            for (iterator i=atomBegin(), e=atomEnd(); i!=e; ++i) {
                atom_t& a = atomFAST(*i);
                a.x = *getter++;
                a.y = *getter++;
                a.z = *getter++;
            }
        }

        /* iterators over element ids */
        iterator atomBegin() const { 
            return iterator(0, _deadatoms.empty() ? NULL : &_deadatoms); 
        }
        iterator atomEnd() const { return iterator(maxAtomId(), NULL); }
        iterator bondBegin() const { 
            return iterator(0, _deadbonds.empty() ? NULL : &_deadbonds); 
        }
        iterator bondEnd() const { return iterator(maxBondId(), NULL); }
        iterator residueBegin() const { 
            return iterator(0, _deadresidues.empty() ? NULL : &_deadresidues); 
        }
        iterator residueEnd() const { return iterator(maxResidueId(), NULL); }
        iterator chainBegin() const { 
            return iterator(0, _deadchains.empty() ? NULL : &_deadchains); 
        }
        iterator chainEnd() const { return iterator(maxChainId(), NULL); }

        /* add an element */
        Id addAtom(Id residue);
        Id addBond(Id i, Id j);
        Id addResidue(Id chain);

        /* If ct is not supplied, add chain to first ct, creating a new one
         * if necessary */
        Id addChain(Id ct=BadId);
        Id addCt();

        /* delete an element */
        void delAtom(Id id);
        void delBond(Id id);
        void delResidue(Id id);
        void delChain(Id id);
        void delCt(Id id);

        /* One more than highest valid id */
        Id maxAtomId() const { return _atoms.size(); }
        Id maxBondId() const { return _bonds.size(); }
        Id maxResidueId() const { return _residues.size(); }
        Id maxChainId() const { return _chains.size(); }
        Id maxCtId() const { return _cts.size(); }

        /* list of elements ids */
        IdList atoms() const;
        IdList bonds() const;
        IdList residues() const;
        IdList chains() const;
        IdList cts() const;

        /* number of elements */
        Id atomCount() const { return _atoms.size() - _deadatoms.size(); }
        Id bondCount() const { return _bonds.size() - _deadbonds.size(); }
        Id residueCount() const {return _residues.size()-_deadresidues.size();}
        Id chainCount() const { return _chains.size() - _deadchains.size(); }
        Id ctCount() const { return _cts.size() - _deadcts.size(); }

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
        Id chainCountForCt(Id id) const {
            if (id>=_ctchains.size()) return 0;
            return _ctchains[id].size();
        }

        /* list of subelements */
        IdList const& bondsForAtom(Id id) const {
            if (id>=_bondindex.size()) return _empty;
            return _bondindex[id];
        }
        /* filtered list of subelements.  Predicate implements
         * bool operator()(bond_t const& b) const; */
        template <typename T>
        IdList filteredBondsForAtom(Id id, T const& predicate) const {
            IdList const& src = _bondindex.at(id);
            IdList dst;
            for (IdList::const_iterator it=src.begin(); it!=src.end(); ++it) {
                Id bid = *it;
                bond_t const& b = bondFAST(bid);
                if (predicate(b)) {
                    dst.push_back(bid);
                }
            }
            return dst;
        }

        IdList const& atomsForResidue(Id id) const {
            if (id>=_residueatoms.size()) return _empty;
            return _residueatoms[id];
        }
        IdList const& residuesForChain(Id id) const {
            if (id>=_chainresidues.size()) return _empty;
            return _chainresidues[id];
        }
        IdList const& chainsForCt(Id id) const {
            if (id>=_ctchains.size()) return _empty;
            return _ctchains[id];
        }

        IdList atomsForCt(Id id) const;
        Id atomCountForCt(Id id) const;
        IdList bondsForCt(Id id) const;
        Id bondCountForCt(Id id) const;

        IdList residuesForCt(Id id) const;


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
        bool hasCt(Id id) const {
            return id<_cts.size() && !_deadcts.count(id);
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
        template <typename T> void delCts(const T& ids) {
            typedef typename T::const_iterator iterator;
            for (iterator i=ids.begin(), e=ids.end(); i!=e; ++i) delCt(*i);
        }

        /* operations on term tables */
        std::vector<String> tableNames() const;
        /* fetch the table with the given name; return NULL if not present */
        TermTablePtr table(const String& name) const;
        /* get the name of the table; throw if table doesn't belong to this */
        String tableName(boost::shared_ptr<TermTable const> table) const;
        /* rename the table with the given name; throw if no such table,
         * or if a table named newname already exists.  */
        void renameTable(String const& oldname, String const& newname);

        TermTablePtr addTable(const String& name, Id natoms,
                              ParamTablePtr ptr = ParamTablePtr() );
        void delTable(const String& name);
        void removeTable(TermTablePtr terms);

        /* invoke coalesce on each table */
        void coalesceTables();

        /* operations on auxiliary tables */
        std::vector<String> auxTableNames() const;
        ParamTablePtr auxTable(String const& name) const;
        void addAuxTable(String const& name, ParamTablePtr aux);
        void delAuxTable(String const& name);
        void removeAuxTable(ParamTablePtr aux);

        /* assign the atom to the given residue */
        void setResidue(Id atom, Id residue);

        /* assign the residue to the given chain */
        void setChain(Id residue, Id chain);

        /* assign the chain to the given ct */
        void setCt(Id chain, Id ct);

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

        /* bonded atoms satisfying a predicate.  predicate implements
         * bool operator()(atom_t const& atm) const; */
        template <typename T>
        IdList filteredBondedAtoms(Id id, T const& predicate) const {
            if (id>=_bondindex.size()) return _empty;
            IdList const& src = _bondindex[id];
            IdList dst;
            for (IdList::const_iterator it=src.begin(); it!=src.end(); ++it) {
                Id bid = *it;
                Id other = bond(bid).other(id);
                atom_t const& atm = atom(other);
                if (predicate(atm)) {
                    dst.push_back(other);
                }
            }
            return dst;
        }
    
        /* Assign atom and residue types; do this after loading a new
         * system from a file or creating it from scratch.  This method
         * also calls updateFragids() for you. */
        void analyze();

        /* update the fragid of each atom according to its bond topology:
         * bonded atoms share the same fragid.  Return the number of
         * frags found, and atomid to fragment partitioning if requested */
        Id updateFragids(MultiIdList* fragments=NULL);

        /* Return ids of atoms based on their order of appearance in
         * a depth-first traversal of the structure hierarchy. */
        IdList orderedIds() const;

        /* Return a pfx graph for the system topology.  Fails if the system
         * has deleted atoms. */
        pfx::Graph* topology() const;
    };

    typedef boost::shared_ptr<System> SystemPtr;

    /* A helper object for importers that manages the logic of 
     * mapping atoms to residues and residue to chains */
    class SystemImporter {
        SystemPtr sys;

        struct ChnKey {
            Id     ct;
            String name;
            String segid;

            ChnKey() {}
            ChnKey(Id c, String const& nm, String const& seg)
            : ct(c), name(nm), segid(seg) {}

            bool operator<(ChnKey const& c) const {
                if (ct!=c.ct) return ct<c.ct;
                int rc = name.compare(c.name);
                if (rc) return rc<0;
                return segid.compare(c.segid)<0;
            }
        };
                
        struct ResKey {
            Id      chain;
            int     resnum;
            String  resname;
            String  insertion;

            ResKey() {}
            ResKey(Id chn, int num, String const& name, String const& insert) 
            : chain(chn), resnum(num), resname(name), insertion(insert) {}

            bool operator<(const ResKey& r) const {
                if (chain!=r.chain) return chain<r.chain;
                if (resnum!=r.resnum) return resnum<r.resnum;
                if (resname!=r.resname) return resname<r.resname;
                return insertion.compare(r.insertion)<0;
            }
        };

        typedef std::map<ChnKey,Id> ChnMap;
        typedef std::map<ResKey,Id> ResMap;
        ResMap resmap;
        ChnMap chnmap;

        Id chnid;
        Id resid;

    public:
        explicit SystemImporter(SystemPtr s) 
        : sys(s), chnid(BadId), resid(BadId) {}

        /* process existing atoms in the system */
        void initialize(IdList const& atoms);

        /* mark a chain as terminated, so that subsequent atoms with
         * the same chain name will be added to a new chain object. */
        bool terminateChain(std::string chain, std::string segid, Id ct=0);

        /* add an atom, after first constructing necessary parent
         * chain and/or residue object.  All string inputs will
         * have leading and trailing whitespace removed. */
        Id addAtom(std::string chain, std::string segid, 
                   int resnum, std::string resname, 
                   std::string atomname,
                   std::string insertion="",
                   Id ct=0);
    };

}}

#endif
