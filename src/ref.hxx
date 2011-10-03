#ifndef desres_msys_ref_hxx
#define desres_msys_ref_hxx

#include "system.hxx"
#include <stdexcept>
#include <sstream>
#include <boost/weak_ptr.hpp>
#include <boost/shared_ptr.hpp>

namespace desres { namespace msys {

    /* Ref classes are intended to operate as the public interface for
     * the pieces that make up a System. */
    template <class Obj, class Parent = System>
    class Ref {
        boost::weak_ptr<Parent> _sys;
        Id      _id;
    
    protected:
        Ref(boost::shared_ptr<Parent> sys, Id id) 
        : _sys(sys), _id(id) {}

    public:
        typedef boost::shared_ptr<Parent> parent_ptr_type;
        typedef std::vector<Obj> list_type;

        /* construct a vector of Obj from the parent and ids */
        template <typename Container>
        static list_type list(parent_ptr_type parent, const Container& ids) {
            list_type v;
            typedef typename Container::const_iterator const_iterator;
            for (const_iterator i=ids.begin(), e=ids.end(); 
                    i!=e; ++i) {
                v.push_back(Obj(parent, *i));
            }
            return v;
        }

        /* return the ids of the objects in the list.  Throws an exception
         * if the parent of each object is not equal to the given parent */
        static IdList ids(const list_type& v, const parent_ptr_type& p) {
            IdList L(v.size());
            for (unsigned i=0, n=v.size(); i<n; i++) {
                L[i]=v[i].id();
                if (v[i].sys() != p) {
                    throw std::runtime_error("ids not from specified parent");
                }
            }
            return L;
        }

        static parent_ptr_type extract(const list_type& v, IdList& ids) {
            ids.resize(v.size());
            if (!v.size()) return parent_ptr_type();
            parent_ptr_type p = v[0].sys();
            for (unsigned i=0, n=v.size(); i<n; i++) {
                ids[i] = v[i].id();
                if (v[i].sys() != p) {
                    throw std::runtime_error("not all elements have same parent");
                }
            }
            return p;
        }

        /* remove this instance of the object */
        void destroy();

        /* does the object still exist? */
        bool exists() const;

        /* name of this object type */
        std::string objtype() const;

        /* id within the structure */
        Id id() const { return _id; }

        /* accessors for structure */
        parent_ptr_type sys() const { 
            if (_sys.expired()) {
                std::stringstream ss;
                ss << "Parent of " << objtype() << " id " << id() 
                   << " no longer exists";
                throw std::runtime_error(ss.str());
            }
            return _sys.lock(); 
        }

        /* ordering comparison for std::map, std::set, etc.  Yes, this could
         * result in non-deterministic behavior if you are depending on the
         * ordering of objects from different systems.  */
        bool operator<(const Obj& obj) const {
            if (id()!=obj.id()) return id()<obj.id();
            return _sys < obj._sys;
        }

        /* identity comparison for python containers.  Here we can use 
         * the structure in the comparison since we are checking for equality. 
         */
        bool operator==(const Obj& obj) const {
            return id()==obj.id() && sys()==obj.sys();
        }
        bool operator!=(const Obj& obj) const {
            return !(*this==obj);
        }
    };

}}

#endif
