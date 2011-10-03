#ifndef desres_msys_view_hxx
#define desres_msys_view_hxx

#include "atom.hxx"
#include "bond.hxx"
#include "residue.hxx"
#include "chain.hxx"

namespace desres { namespace msys {

    class AtomView;
    class BondView;
    class ResidueView;
    class ChainView;
    class View;

    typedef std::vector<AtomView>     AtomViewList;
    typedef std::vector<BondView>     BondViewList;
    typedef std::vector<ResidueView>  ResidueViewList;
    typedef std::vector<ChainView>    ChainViewList;

    class View {
        System* _sys;
    
        typedef std::set<Id> IdSet;
        IdSet   _atoms;
        IdSet   _bonds;
        IdSet   _residues;
        IdSet   _chains;
    
    
    public:
        View( System* system, const IdList& atoms );
        System* sys() { return _sys; }
    
        Id atomCount() const { return _atoms.size(); }
        AtomViewList atoms() { return AtomViewList(this, _atoms); }

        AtomViewList atoms(ResidueView res);
        Id atomCount(ResidueView res) const;

        Id bondCount() const { return _bonds.size(); }
        BondViewList bonds() { return BondViewList(this, _bonds); }

        BondViewList bonds(AtomView atom);
        Id atomCount(AtomView atom) const;

    };

    class AtomView : public Atom {
        View* _view;
        /* could cache things like bond list, residue id, etc. */
        
    public:
        AtomView(View* view, Id id) 
        : Atom(view->sys(), id), _view(view) {}

        Id bondCount() const {
            return _view->atomCount(*this);
        }
        BondViewList bonds() {
            return _view->bonds(*this);
        }
    };

}}

#endif
