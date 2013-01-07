#ifndef desres_viparr_bondFilters_hxx
#define desres_viparr_bondFilters_hxx

#include <msys/system.hxx>

namespace desres { namespace msys {
   
    class bondFilter {
        SystemPtr _sys;
    public:
        bondFilter(SystemPtr sys) : _sys(sys) { }
        virtual ~bondFilter(){};
        /* Keep bonded atom of a bond if valid */
        bool operator()(const msys::atom_t& a) const{
           return atomValid(a);
        }

        /* Keep a bond if both atoms should be kept */
        bool operator()(const msys::bond_t& b) const {
            return atomValid(_sys->atom(b.i)) && atomValid(_sys->atom(b.j));
        }
  
        /* return true if atom should be kept */ 
        virtual bool atomValid(const msys::atom_t& a)const=0;
    };


    /* Exclude bonded atoms that are virtuals and bonds to virtuals */
    class bondedVirtualsFilter: public bondFilter {
    public:
        bondedVirtualsFilter(SystemPtr sys) : bondFilter(sys) { }
        virtual bool atomValid(const msys::atom_t& a) const {
            return (a.atomic_number > 0);
        }

    };

    /* Exclude bonded atoms that are virtuals/metals and bonds to virtuals/metals */
    class bondedVirtualsAndMetalsFilter : public bondFilter{
    public:
        bondedVirtualsAndMetalsFilter(SystemPtr sys) : bondFilter(sys) { }
        virtual bool atomValid(const msys::atom_t& a) const {
            bool virtualOrMetal=(
                                  (a.atomic_number < 1) ||
                                  (a.atomic_number == 13) ||
                                  (a.atomic_number >= 21 && a.atomic_number <= 32) ||
                                  (a.atomic_number >= 39 && a.atomic_number <= 51) ||
                                  (a.atomic_number >= 57 && a.atomic_number <= 84) ||
                                  (a.atomic_number >= 89 && a.atomic_number <= 117)
                                  );
            return !virtualOrMetal;
        }
    };

}}

#endif
