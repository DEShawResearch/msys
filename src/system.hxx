#ifndef desres_msys_system_hxx
#define desres_msys_system_hxx

#include <vector>
#include <map>
#include <boost/enable_shared_from_this.hpp>

#include "provenance.hxx"
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

        explicit Vec3(const Float* xyz)
        : x(xyz[0]), y(xyz[1]), z(xyz[2]) {}

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

        Vec3& operator-=(Vec3 const& o) {
            x-=o.x;
            y-=o.y;
            z-=o.z;
            return *this;
        }

        Vec3& operator *=(Float s) {
            x *= s;
            y *= s;
            z *= s;
            return *this;
        }

        Float dot(Vec3 const& o) const {
            return x*o.x + y*o.y + z*o.z;
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
        String es_funct;
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
    
        String name;
        AtomType type;
        Float resonant_charge;
    
        atom_t() 
        : fragid(BadId), residue(BadId), atomic_number(0), formal_charge(0),
          x(0), y(0), z(0), charge(0), vx(0), vy(0), vz(0), mass(0), type(),
          resonant_charge()
        {}
    };
    
    struct bond_t {
        Id  i;                  /* id of first bond partner     */
        Id  j;                  /* id of second bond partner    */
        int order;              /* formal bond order            */
        Float resonant_order;   /* resonant bond order          */
    
        bond_t() : i(BadId), j(BadId), order(1), resonant_order(1) {}
        bond_t(Id ai, Id aj) : i(ai), j(aj), order(1), resonant_order(1) {}
        Id other(Id id) const { return id==i ? j : i; }
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
        String  name;
        String insertion;
        ResidueType type;
    
        residue_t() : chain(BadId), resid(), type() {}
    };
    
    struct chain_t {
        String  name;
        String  segid;
    
        chain_t() {}
    };

    typedef std::pair<Id,Id> glue_t;

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

        /* auxiliary tables.  basically a hack for cmap */
        typedef std::map<String, ParamTablePtr> AuxTableMap;
        AuxTableMap _auxtables;

        /* atom selection macros */
        typedef std::map<String, String> MacroMap;
        MacroMap    _macros;

        /* the glue table */
        typedef std::set<glue_t> GlueSet;
        GlueSet     _glue;

        /* provenance.  Ideally, you would append to this just before 
         * serializing to disk. */
        std::vector<Provenance> _provenance;

        /* cached SSSR */
        boost::shared_ptr<MultiIdList> _allRelevantSSSR;

        /* create only as shared pointer. */
        System();
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
        /* filtered list of subelements.  Predicate implements
         * bool operator()(bond_t const& b) const; */
        template <typename T>
        IdList filteredBondsForAtom(Id id, T const& predicate) const {
            IdList const& src = _bondindex[id];
            IdList dst;
            for (IdList::const_iterator it=src.begin(); it!=src.end(); ++it) {
                Id bid = *it;
                bond_t const& b = bond(bid);
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

        /* return the cached result of a call to GetSSSR(this, atoms(), true).
         * The cache will be invalidated only by a call to updateFragids(). */
        boost::shared_ptr<MultiIdList> allRelevantSSSR();

        /* Return ids of atoms based on their order of appearance in
         * a depth-first traversal of the structure hierarchy. */
        IdList orderedIds() const;

        /* atom selection macros are strings that are expanded at parse time
         * into longer strings.  */
        void addSelectionMacro(std::string const& macro,
                               std::string const& definition);

        /* Remove the given macro */
        void delSelectionMacro(std::string const& macro);

        /* Get the definition corresponding to the given macro.  Returns empty
         * string if not found. */
        std::string const& selectionMacroDefinition(std::string const& m) const;

        /* number of defined selection macros */
        Id selectionMacroCount() const;

        /* Return a list of all the selection macros */
        std::vector<std::string> selectionMacros() const;

        /* Remove all macros */
        void clearSelectionMacros();

        /* Assign the default selection macros */
        void initSelectionMacros();

        /* copy macros from another system */
        void copySelectionMacros(System const& m);


        /**** glue ***/

        /* number of glue pairs */
        Id glueCount() const;

        /* the glue pairs */
        std::vector<glue_t> gluePairs() const;

        /* glue present? */
        bool hasGluePair(Id p0, Id p1) const;

        /* add a glue pair.  Return true if a pair was added; false if
         * the pair was already present */
        bool addGluePair(Id p0, Id p1);

        /* remove a glue pair.  Return true if found, false if not. */
        bool delGluePair(Id p0, Id p1);
    };

    typedef boost::shared_ptr<System> SystemPtr;

    /* A helper object for importers that manages the logic of 
     * mapping atoms to residues and residue to chains */
    class SystemImporter {
        SystemPtr sys;

        struct ChnKey {
            String name;
            String segid;

            ChnKey() {}
            ChnKey(String const& nm, String const& seg)
            : name(nm), segid(seg) {}

            bool operator<(ChnKey const& c) const {
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

        /* add an atom, after first constructing necessary parent
         * chain and/or residue object.  All string inputs will
         * have leading and trailing whitespace removed. */
        Id addAtom(std::string chain, std::string segid, 
                   int resnum, std::string resname, 
                   std::string atomname,
                   std::string insertion="");
    };

}}

#endif
