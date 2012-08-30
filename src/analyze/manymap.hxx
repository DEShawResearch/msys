#ifndef MANYMAP_H
#define MANYMAP_H 1

#include <map>
#include <vector>


/* Convience wrapper around a map of vectors (as an alternative to std::multimap). 
   Benefits over std::multimap and/or std::map of std::vectors include:
   1) size reports number of unique keys (not total number of elements as in multimap)
   2) count reports number of values with given key (not 1 as in map of vectors)
   3) To iterate over the set of unique keys, std::multimap requires you to either keep 
   track of them (to use equal range) or keep track of the last seen key (when using iterators)
   which is just messy. Here, every change in the iterator is a new key... Very clean.
   4) equal_range returns the iterators to the underlying vector of values
   5) Ive noticed pretty substantial performance increases when using map of vectors versus 
   multimap.

   supports some c++11 features

   TODO: move constructors and move initalizers
         list initalizers
         emplace, emplace hint
*/

namespace desres {

    template <typename Key, typename T, typename Compare = std::less<Key>,
              typename Allocator = std::allocator<std::pair<const Key, std::vector<T> > > >
    class manymap {
    public:
        typedef Key                                       key_type;
        typedef T                                         entry_type;
        typedef std::vector<T>                            mapped_type;
        typedef std::pair<const Key, std::vector<T> >     value_type;
        typedef Compare                                   key_compare;
        typedef Allocator                                 allocator_type;
        typedef std::pair<const Key,  T >                 insert_type;

    private:
        typedef std::map< key_type, mapped_type, key_compare, allocator_type > manymap_t;
        manymap_t _manymap;
 
    public:
        typedef typename manymap_t::value_compare          value_compare;
        typedef typename manymap_t::reference              reference;
        typedef typename manymap_t::const_reference        const_reference;
        typedef typename manymap_t::iterator               iterator;
        typedef typename manymap_t::const_iterator         const_iterator;
        typedef typename manymap_t::size_type              size_type;
        typedef typename manymap_t::difference_type        difference_type;
        typedef typename manymap_t::pointer                pointer;
        typedef typename manymap_t::const_pointer          const_pointer;
        typedef typename manymap_t::reverse_iterator       reverse_iterator;
        typedef typename manymap_t::const_reverse_iterator const_reverse_iterator;

        // Constructors
        explicit manymap(const Compare& comp = Compare(), 
                         const Allocator& alloc = Allocator())
            : _manymap(comp, alloc) { }

        template<typename InputIterator>
        manymap(InputIterator first, InputIterator last,
                const Compare& comp = Compare(),
                const Allocator& alloc = Allocator())
            : _manymap(comp, alloc) { 
            insert(first, last); 
        }
        manymap(const manymap_t& x) : _manymap(x._manymap) {}

        manymap_t& operator=(const manymap_t& x) { 
            _manymap = x._manymap;
            return *this;
        }

        //  Allocator:
        inline allocator_type get_allocator() const { 
            return _manymap.get_allocator(); 
        }


        // Iterators:
        inline iterator begin(){ 
            return _manymap.begin(); 
        }
        inline const_iterator begin() const { 
            return _manymap.begin(); 
        }
        inline iterator end(){ 
            return _manymap.end(); 
        }
        inline const_iterator end() const { 
            return _manymap.end(); 
        }
        inline reverse_iterator rbegin(){ 
            return _manymap.rbegin(); 
        }
        inline const_reverse_iterator rbegin() const { 
            return _manymap.rbegin(); 
        }
        inline reverse_iterator rend(){ 
            return _manymap.rend(); 
        }
        inline const_reverse_iterator rend() const { 
            return _manymap.rend(); 
        }
        
        inline const_iterator cbegin() const { 
            return _manymap.begin(); 
        }
        inline const_iterator cend() const { 
            return _manymap.end(); 
        }
        inline const_reverse_iterator crbegin() const { 
            return _manymap.rbegin(); 
        }
        inline const_reverse_iterator crend() const { 
            return _manymap.rend(); 
        }

        // Capacity:
        inline bool empty() const { 
            return _manymap.empty(); 
        }
        inline size_type size() const { 
            return _manymap.size(); 
        }
        inline size_type max_size() const { 
            return _manymap.max_size(); 
        }

        // Element access:
        inline mapped_type& operator[](const key_type& x) { 
            return _manymap[x]; 
        }
        inline mapped_type& at(const key_type& x){
            return _manymap.at(x);
        }
        const mapped_type& at(const key_type& x) const {
            return _manymap.at(x);
        }

        // Modifiers:
        /* Single insert_type insert */
        std::pair<iterator, bool> insert(const insert_type& x){
            iterator iter=_manymap.lower_bound(x.first);
            if(iter ==_manymap.end() || _manymap.key_comp()(x.first, 
                                                            iter->first)){
                iter=_manymap.insert(iter, typename manymap_t::value_type(x.first, mapped_type()) );
            }
            iter->second.push_back(x.second);
            return std::make_pair(iter,true);
        }

        /* Iterator insert */
        template<typename InputIterator>
        void insert(InputIterator first, InputIterator last){
            for(;first != last; ++first){
                iterator iter=_manymap.lower_bound(first->first);
                if(iter == _manymap.end() || _manymap.key_comp()(first->first,
                                                                 iter->first)){
                    iter=_manymap.insert(iter,*first);
                }else{
                    iter->second.insert(iter->second.end(),
                                        first->second.begin(),first->second.end());
                }
            }
        }

        inline void erase(iterator position){
            _manymap.erase(position);
        }

        inline iterator erase(const_iterator position){
            iterator result(position);
            ++result;
            _manymap.erase(position);
            return result;
        }

        inline size_type erase(const key_type& x) {
            iterator position=_manymap.find(x);
            if(position==_manymap.end()) return 0;
            size_type result=position->second.size();
            _manymap.erase(position);
            return result; 
        }
        inline void erase(iterator first, iterator last) { 
            _manymap.erase(first, last); 
        }

        inline iterator erase(const_iterator first, const_iterator last) { 
            iterator result(last);
            _manymap.erase(first,last);
            return result;            
        }

        inline void swap(manymap_t& x) { 
            _manymap.swap(x._manymap); 
        }
        inline void clear() { 
            _manymap.clear(); 
        }

        // Observers:
        inline key_compare key_comp() const { 
            return _manymap.key_comp(); 
        }
        inline value_compare value_comp() const { 
            return _manymap.value_comp(); 
        }

        // Operations:
        inline iterator find(const key_type& x) { 
            return _manymap.find(x); 
        }
        inline const_iterator find(const key_type& x) const { 
            return _manymap.find(x); 
        }
        size_type count(const key_type& x) const { 
            iterator iter=_manymap.find(x);
            if(iter==_manymap.end()) return 0;
            return iter->second.size();
        }

        inline iterator lower_bound(const key_type& x) { 
            return _manymap.lower_bound(x); 
        }
        inline const_iterator lower_bound(const key_type& x) const { 
            return _manymap.lower_bound(x); 
        }
        inline iterator upper_bound(const key_type& x) { 
            return _manymap.upper_bound(x); 
        }
        inline const_iterator upper_bound(const key_type& x) const { 
            return _manymap.upper_bound(x); 
        }
        inline std::pair<typename mapped_type::iterator, 
                         typename mapped_type::iterator>
        equal_range(const key_type& x) { 
            iterator iter=_manymap.find(x);
            if(iter==_manymap.end()){
              typename mapped_type::iterator invalid(mapped_type::pointer(0));
              return std::pair<typename mapped_type::iterator, 
                               typename mapped_type::iterator>(invalid,
                                                               invalid);
            } 
            return std::pair<typename mapped_type::iterator, 
                             typename mapped_type::iterator>(iter->second.begin(),
                                                             iter->second.end());
        }
        inline std::pair<typename mapped_type::const_iterator,
                         typename mapped_type::const_iterator>
        equal_range(const key_type& x) const {
            iterator iter=_manymap.find(x);
            if(iter==_manymap.end()){
              typename mapped_type::iterator invalid(mapped_type::pointer(0));
              return std::pair<typename mapped_type::const_iterator,
                               typename mapped_type::const_iterator>(invalid,
                                                                    invalid);
            }
            return std::pair<typename mapped_type::const_iterator,
                             typename mapped_type::const_iterator>(iter->second.begin(),
                                                                   iter->second.end());
        }

        template <typename _Key, typename _T, typename _Compare, typename _Allocator>
        friend bool operator==(const manymap<_Key,_T,_Compare,_Allocator>& x,
                               const manymap<_Key,_T,_Compare,_Allocator>& y);
        template <typename _Key, typename _T, typename _Compare, typename _Allocator>
        friend bool operator< (const manymap<_Key,_T,_Compare,_Allocator>& x,
                               const manymap<_Key,_T,_Compare,_Allocator>& y);


    };

    template <typename Key, typename T, typename Compare, typename Allocator>
    inline bool operator==(const manymap<Key,T,Compare,Allocator>& x,
                    const manymap<Key,T,Compare,Allocator>& y){
        return (x._manymap == y._manymap);
    }
    template <typename Key, typename T, typename Compare, typename Allocator>
    inline bool operator< (const manymap<Key,T,Compare,Allocator>& x,
                           const manymap<Key,T,Compare,Allocator>& y){
        return (x._manymap < y._manymap);
    }
    template <typename Key, typename T, typename Compare, typename Allocator>
    inline bool operator!=(const manymap<Key,T,Compare,Allocator>& x,
                           const manymap<Key,T,Compare,Allocator>& y){
        return !(x == y);
    }
    template <typename Key, typename T, typename Compare, typename Allocator>
    inline bool operator> (const manymap<Key,T,Compare,Allocator>& x,
                           const manymap<Key,T,Compare,Allocator>& y){
        return (y < x);
    }
    template <typename Key, typename T, typename Compare, typename Allocator>
    inline bool operator>=(const manymap<Key,T,Compare,Allocator>& x,
                           const manymap<Key,T,Compare,Allocator>& y){
        return !(x < y);
    }
    template <typename Key, typename T, typename Compare, typename Allocator>
    inline bool operator<=(const manymap<Key,T,Compare,Allocator>& x,
                           const manymap<Key,T,Compare,Allocator>& y){
        return !(y < x);
    }

    template <typename Key, typename T, typename Compare, typename Allocator>
    inline void swap(manymap<Key,T,Compare,Allocator>& x,
                     manymap<Key,T,Compare,Allocator>& y){ 
        x.swap(y);
    }

}

#endif
