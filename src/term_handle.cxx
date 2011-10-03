#include "term_handle.hxx"

using namespace desres::msys;

namespace desres { namespace msys {
    template<> std::string Ref<Term,TermTable>::objtype() const { 
        return "Term"; 
    }
    template <> bool Ref<Term,TermTable>::exists() const { 
        return sys()->hasTerm(id()); 
    }
}}

Term TermTableHandle::addTerm(const AtomList& atoms, Param param ) {
    Id id=param.id();
    if (id!=BadId && param.params()!=params()) {
        throw std::runtime_error(
                "term table uses different param table than param");
    }
    return Term(_ptr, _ptr->addTerm(Atom::ids(atoms, _ptr->system()), id));
}

TermList TermTableHandle::terms() {
    return Term::list( _ptr, _ptr->terms());
}
