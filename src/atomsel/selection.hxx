#ifndef desres_msys_atomsel_selection_hxx
#define desres_msys_atomsel_selection_hxx

#include "../types.hxx"
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

namespace desres { namespace msys { namespace atomsel {

    /* A selection is a list of booleans */

    class Selection {

        std::vector<unsigned char> flags;
        Id _size;

        void alloc() { 
            flags.resize(_size + sizeof(size_t)-1);
        }

    public:
        /* initialize with size and list of flags starting out as true */
        Selection(Id size, const IdList& ids) : _size(size) {
            alloc();
            IdList::const_iterator i, e;
            for (i=ids.begin(), e=ids.end(); i!=e; ++i) flags[*i]=1;
        }

        /* construct an empty selection of the given size */
        explicit Selection(Id size) : _size(size) {
            alloc();
        }
        
        /* clear all flags */
        void clear() { memset(&flags[0], 0, _size); }

        /* number selected */
        Id count() const {
            Id cnt=0;
            for (Id i=0, n=_size; i<n; i++) cnt += flags[i];
            return cnt;
        }

        /* total number of flags */
        inline Id size() const { return _size; }

        /* accessors */
        inline unsigned char&       operator[](Id i)       { return flags[i]; }
        inline unsigned char const& operator[](Id i) const { return flags[i]; }

        /* return the selected ids */
        IdList ids() const {
            IdList tmp;
            for (Id i=0, n=_size; i!=n; i++) if (flags[i]) tmp.push_back(i);
            return tmp;
        }

        /* AND a selection into this one */
        Selection& intersect(Selection const& other) {
            size_t* mine = (size_t*)&flags[0];
            const size_t* that = (const size_t*)&other.flags[0];
            Id i,n = flags.size()/sizeof(size_t);
            for (i=0; i<n; i++) mine[i] &= that[i];
            return *this;
        }

        /* OR a selection into this one */
        Selection& add(Selection const& other) {
            size_t* mine = (size_t*)&flags[0];
            const size_t* that = (const size_t*)&other.flags[0];
            Id i,n = flags.size()/sizeof(size_t);
            for (i=0; i<n; i++) mine[i] |= that[i];
            return *this;
        }

        /* AND NOT a selection into this one (subtract) */
        Selection& subtract(Selection const& other) {
            size_t* mine = (size_t*)&flags[0];
            const size_t* that = (const size_t*)&other.flags[0];
            Id i,n = flags.size()/sizeof(size_t);
            for (i=0; i<n; i++) mine[i] &= ~(that[i]);
            return *this;
        }
    };

}}}

#endif
