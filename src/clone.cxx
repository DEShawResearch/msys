#include "clone.hxx"
#include "append.hxx"
#include <algorithm>
#include <stdexcept>
#include <sstream>

using namespace desres::msys;

namespace {
    /* a predicate that returns true if a term is not in a subset */
    struct BadTerm {
        TermTablePtr _table;
        IdList const& _atmmap;

        BadTerm( TermTablePtr table, IdList const& atmmap ) 
        : _table(table), _atmmap(atmmap) {}
        bool operator()( Id const& id ) {
            IdList atoms = _table->atoms(id);
            for (unsigned i=0; i<atoms.size(); i++) {
                if (bad(_atmmap[atoms[i]])) return true;
            }
            return false;
        }
    };
}

SystemPtr desres::msys::Clone( SystemPtr src, IdList const& atoms ) {

    /* check for duplicates */
    {
        IdList tmp(atoms);
        std::sort(tmp.begin(), tmp.end());
        if (std::unique(tmp.begin(), tmp.end())!=tmp.end()) {
            throw std::runtime_error(
                    "Clone: atoms argument contains duplicates");
        }
    }
    SystemPtr dst = System::create();

    /* copy selection macros */
    dst->copySelectionMacros(*src);

    /* Mappings from src ids to dst ids */
    IdList atmmap(src->maxAtomId(), BadId);
    IdList resmap(src->maxResidueId(), BadId);
    IdList chnmap(src->maxChainId(), BadId);

    /* copy atom properties */
    Id nprops = src->atomPropCount();
    IdList propmap(nprops);
    for (Id i=0; i<nprops; i++) {
        propmap[i] = dst->addAtomProp(
                src->atomPropName(i), src->atomPropType(i));
    }

    /* copy bond properties */
    Id nbprops = src->bondPropCount();
    IdList bpropmap(nbprops);
    for (Id i=0; i<nbprops; i++) {
        bpropmap[i] = dst->addBondProp(
                src->bondPropName(i), src->bondPropType(i));
    }

    /* Build structure for subset of atoms */
    for (Id i=0; i<atoms.size(); i++) {
        Id srcatm = atoms[i];
        Id srcres = src->atom(srcatm).residue;
        Id srcchn = src->residue(srcres).chain;

        if (!src->hasAtom(srcatm)) {
            std::stringstream ss;
            ss << "Clone: atoms argument contains deleted atom id " 
               << srcatm;
            throw std::runtime_error(ss.str());
        }

        Id dstchn = chnmap[srcchn];
        Id dstres = resmap[srcres];

        if (bad(dstres)) {
            if (bad(dstchn)) {
                dstchn = chnmap[srcchn] = dst->addChain();
                dst->chain(dstchn) = src->chain(srcchn);
            }
            dstres = resmap[srcres] = dst->addResidue(dstchn);
            dst->residue(dstres) = src->residue(srcres);
            dst->residue(dstres).chain = dstchn;
        }
        Id dstatm = dst->addAtom(dstres);
        atmmap[srcatm] = dstatm;
        /* Copy built-in properties */
        dst->atom(dstatm) = src->atom(srcatm);
        /* Restore the overwritten residue id */
        dst->atom(dstatm).residue = dstres;
        /* Copy additional atom properties */
        for (Id j=0; j<nprops; j++) {
            dst->atomPropValue(dstatm,propmap[j]) = 
            src->atomPropValue(srcatm, j);
        }
    }

    /* Build bonds whose atoms are fully within the subset */
    for (Id i=0; i<atoms.size(); i++) {
        Id srcatm = atoms[i];
        IdList const& bonds = src->bondsForAtom(srcatm);
        for (Id j=0; j<bonds.size(); j++) {
            bond_t const& srcbond = src->bond(bonds[j]);
            Id srcother = srcbond.other(srcatm);
            Id dstother = atmmap[srcother];
            if (!bad(dstother)) {
                Id dstatm = atmmap[srcatm];
                Id dstbnd = dst->addBond(dstatm, dstother);
                dst->bond(dstbnd).order = srcbond.order;

                /* Copy additional bond properties */
                for (Id k=0; k<nbprops; k++) {
                    dst->bondPropValue(dstbnd,bpropmap[k]) = 
                    src->bondPropValue(bonds[j], k);
                }
            }
        }
    }

    /* Add term tables */
    std::vector<String> tablenames = src->tableNames();
    for (unsigned i=0; i<tablenames.size(); i++) {
        std::string const& name = tablenames[i];
        TermTablePtr srctable = src->table(name);
        TermTablePtr dsttable = dst->addTable(name, srctable->atomCount());
        dsttable->category = srctable->category;
        IdList terms = srctable->terms();
        IdList::iterator iter = std::remove_if( 
                terms.begin(), terms.end(), BadTerm(srctable, atmmap));
        terms.erase(iter, terms.end());
        AppendTerms( dsttable, srctable, atmmap, terms );
    }

    /* copy all global data */
    dst->name = src->name;
    dst->global_cell = src->global_cell;
    dst->nonbonded_info = src->nonbonded_info;

    /* add/replace extra tables */
    std::vector<String> extras = src->auxTableNames();
    for (unsigned i=0; i<extras.size(); i++) {
        std::string const& name = extras[i];
        ParamTablePtr srcparams = src->auxTable(name);
        ParamTablePtr dstparams = ParamTable::create();
        AppendParams( dstparams, srcparams, srcparams->params() );
        dst->addAuxTable( name, dstparams );
    }
    dst->updateFragids();
    return dst;
}

