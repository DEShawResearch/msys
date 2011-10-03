#ifndef desres_msys_term_handle_hxx
#define desres_msys_term_handle_hxx

#include "term_table.hxx"
#include "atom.hxx"
#include "param_handle.hxx"
#include "append.hxx"

namespace desres { namespace msys {

    class TermTableHandle;

    class Term : public Ref<Term, TermTable> {
    public:
        Term( TermTablePtr table=TermTablePtr(), Id id=BadId ) 
        : Ref<Term,TermTable>(table,id) {}

        Atom atom(Id i) {
            IdList ids(sys()->atoms(id()));
            return Atom(sys()->system(), ids.at(i));
        }
        AtomList atoms() {
            return Atom::list(sys()->system(), sys()->atoms(id()));
        }
        Param param() {
            return Param(sys()->paramTable().get(), sys()->param(id()));
        }
        void setParam(Param p) {
            sys()->setParam(id(), p.id());
        }
        Param paramB() {
            return Param(sys()->paramTable().get(), sys()->paramB(id()));
        }
        void setParamB(Param p) {
            sys()->setParamB(id(), p.id());
        }
        inline TermTableHandle table();

        ValueRef prop(Id index) {
            return sys()->termPropValue(id(), index);
        }
        ValueRef prop(const String& key) {
            Id col = sys()->termPropIndex(key);
            if (bad(col)) {
                std::stringstream ss;
                ss << "term table has no term property '" << key << "'";
                throw std::runtime_error(ss.str());
            }
            return prop(col);
        }
        bool hasProp(const String& key) const {
            return !bad(sys()->termPropIndex(key));
        }
    };

    typedef std::vector<Term> TermList;

    class TermTableHandle {
        TermTablePtr _ptr;
    public:
        /* default constructor creates invalid handle */
        TermTableHandle() {}

        /* convertible from TermTablePtr */
        TermTableHandle(TermTablePtr ptr)
        : _ptr(ptr) 
        {}

        ///* constructable from TermList */
        //explicit TermTableHandle(TermList& terms) 
        //: _ptr(terms.sys()->shared_from_this())
        //{}

        /* the underlying shared pointer */
        TermTablePtr ptr() { return _ptr; }

        /* set/get category */
        std::string category() const { return _ptr->category; }
        void setCategory(std::string const& cat) { _ptr->category=cat; }

        /* parameter table */
        ParamTableHandle params() { return _ptr->paramTable(); }

        /* term properties */
        Id addTermProp( const String& name, ValueType type) {
            return _ptr->addTermProp(name,type);
        }
        Id termPropCount() const {
            return _ptr->termPropCount();
        }
        String termPropName(Id i) const {
            return _ptr->termPropName(i);
        }
        ValueType termPropType(Id i) const { 
            return _ptr->termPropType(i); 
        }
        Id termPropIndex(const String& name) {
            return _ptr->termPropIndex(name);
        }

        /* forward TermTable methods */
        Id atomCount() const { return _ptr->atomCount(); }
        Id termCount() const { return _ptr->termCount(); }

        bool operator==(const TermTableHandle& other) const {
            return _ptr==other._ptr;
        }

        /* all the terms */
        TermList terms();

        /* add a term using an atom list and a Param */
        Term addTerm(const AtomList& atoms, Param param);

        IdList appendTerms( TermTableHandle src,
                            IdList const& idmap,
                            IdList const& terms ) {
            return AppendTerms( ptr(), src.ptr(), idmap, terms );
        }
    };

    TermTableHandle Term::table() {
        return sys()->shared_from_this();
    }


}}

#endif

