#ifndef desres_msys_atomsel_selection_hxx
#define desres_msys_atomsel_selection_hxx

#include "../types.hxx"

namespace desres { namespace msys { namespace atomsel {

    /* A selection is a list of booleans */

    class Selection {

        typedef Id flag_type;
        typedef IdList flag_list;
        flag_list flags;

    public:
        /* initialize with size and list of flags starting out as true */
        Selection(Id size, const IdList& ids);

        /* construct an empty selection of the given size */
        Selection(Id size) : flags(size) {}

        /* clear all flags */
        void clear();

        /* number selected */
        Id count() const;

        /* total number of flags */
        inline Id size() const { return flags.size(); }

        /* accessors */
        inline Id&        operator[](Id i)       { return flags[i]; }
        inline Id  const& operator[](Id i) const { return flags[i]; }

        /* return the selected ids */
        IdList ids() const;

        /* AND a selection into this one */
        Selection& intersect(const Selection& other);

        /* OR a selection into this one */
        Selection& add(const Selection& other);

        /* AND NOT a selection into this one (subtract) */
        Selection& subtract(const Selection& other);
    };

}}}

#endif
