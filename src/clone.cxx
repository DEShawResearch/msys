#include "clone.hxx"
#include "append.hxx"
#include <algorithm>
#include <stdexcept>
#include <sstream>
#include <unordered_map>

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

SystemPtr desres::msys::Clone( SystemPtr src, IdList const& atoms,
                               CloneOption::Flags flags) {

    SystemPtr dst = System::create();

    /* Mappings from src ids to dst ids */
    IdList atmmap(src->maxAtomId(), BadId);
    IdList resmap(src->maxResidueId(), BadId);
    IdList chnmap(src->maxChainId(), BadId);
    IdList ctmap(src->maxCtId(), BadId);

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
    for (Id srcatm : atoms) {

        if (!src->hasAtom(srcatm)) {
            MSYS_FAIL("atoms argument contains deleted atom id " << srcatm);
        }
        if (atmmap[srcatm] != BadId) {
            MSYS_FAIL("duplicate atom id " << srcatm);
        }
        auto& atm = src->atomFAST(srcatm);
        if ((flags & CloneOption::StructureOnly) && atm.atomic_number < 1) {
            continue;
        }

        Id srcres = atm.residue;
        Id srcchn = src->residue(srcres).chain;
        Id srcct  = src->chain(srcchn).ct;

        Id dstct  = ctmap[srcct];
        Id dstchn = chnmap[srcchn];
        Id dstres = resmap[srcres];

        if (bad(dstres)) {
            if (bad(dstchn)) {
                if (bad(dstct)) {
                    dstct = ctmap[srcct] = dst->addCt();
                    dst->ct(dstct) = src->ct(srcct);
                }
                dstchn = chnmap[srcchn] = dst->addChain(dstct);
                dst->chain(dstchn) = src->chain(srcchn);
                dst->chain(dstchn).ct = dstct;
            }
            dstres = resmap[srcres] = dst->addResidue(dstchn);
            dst->residue(dstres) = src->residue(srcres);
            dst->residue(dstres).chain = dstchn;
        }
        Id dstatm = dst->addAtom(dstres);
        atmmap[srcatm] = dstatm;
        /* Copy built-in properties */
        dst->atom(dstatm) = atm;
        /* Restore the overwritten residue id */
        dst->atom(dstatm).residue = dstres;
        /* Copy additional atom properties */
        for (Id j=0; j<nprops; j++) {
            dst->atomPropValue(dstatm,propmap[j]) = 
            src->atomPropValue(srcatm, j);
        }
    }

    /* Build bonds whose atoms are fully within the subset */
    for (Id srcatm=0, n=atmmap.size(); srcatm<n; srcatm++) {
        if (bad(atmmap[srcatm])) continue;
        IdList const& bonds = src->bondsForAtom(srcatm);
        for (Id j=0; j<bonds.size(); j++) {
            bond_t const& srcbond = src->bond(bonds[j]);
            Id srcother = srcbond.other(srcatm);
            Id dstother = atmmap[srcother];
            if (!bad(dstother)) {
                Id dstatm = atmmap[srcatm];
                Id dstbnd = dst->addBond(dstatm, dstother);
                dst->bond(dstbnd).order = srcbond.order;
                dst->bond(dstbnd).stereo = srcbond.stereo;
                dst->bond(dstbnd).aromatic = srcbond.aromatic;

                /* Copy additional bond properties */
                for (Id k=0; k<nbprops; k++) {
                    dst->bondPropValue(dstbnd,bpropmap[k]) = 
                    src->bondPropValue(bonds[j], k);
                }
            }
        }
    }

    if (flags & CloneOption::StructureOnly) {
        // skip term tables
        //
    } else if (flags & CloneOption::ShareParams) {
        TermTablePtr srctable, dsttable;
        for (String name : src->tableNames()) {
            srctable = src->table(name);
            dsttable = dst->addTable(name, 
                                     srctable->atomCount(),
                                     srctable->params());
            dsttable->category = srctable->category;
            IdList terms = srctable->terms();
            IdList::iterator iter = std::remove_if(
                    terms.begin(), terms.end(), BadTerm(srctable, atmmap));
            terms.erase(iter, terms.end());
            AppendTerms( dsttable, srctable, atmmap, terms );
        }
        
    } else {

        /* Detect when term tables share a ParamTable.  First pass: sort
         * by param table.  */
        typedef std::vector<std::string> StringList;
        typedef std::unordered_map<ParamTablePtr, StringList> TableMap;
        TableMap tables;
        TableMap::const_iterator it;
        std::vector<String> tablenames = src->tableNames();
        for (unsigned i=0; i<tablenames.size(); i++) {
            std::string const& name = tablenames[i];
            TermTablePtr srctable = src->table(name);
            tables[srctable->params()].push_back(name);
        }

        /* process the unshared tables. */
        for (it=tables.begin(); it!=tables.end(); ++it) {
            if (it->second.size()>1) continue;
            std::string const& name = it->second.at(0);
            TermTablePtr srctable = src->table(name);
            TermTablePtr dsttable = dst->addTable(name, srctable->atomCount());
            dsttable->category = srctable->category;
            IdList terms;
            if (flags & CloneOption::UseIndex) {
                terms = srctable->findWithOnly(atoms);
            } else {
                terms = srctable->terms();
                auto iter = std::remove_if(
                        terms.begin(), terms.end(), BadTerm(srctable, atmmap));
                terms.erase(iter, terms.end());
            }
            AppendTerms( dsttable, srctable, atmmap, terms );
        }

        /* process the tables with shared params */
        for (it=tables.begin(); it!=tables.end(); ++it) {
            StringList const& s = it->second;
            if (it->second.size()<2) continue;
            ParamTablePtr dstparams = ParamTable::create();
            ParamTablePtr srcparams = src->table(s[0])->params();
            IdList src2dst(srcparams->paramCount(), BadId);
            for (unsigned i=0; i<s.size(); i++) {
                std::string const& name = s[i];
                TermTablePtr srctable = src->table(name);
                TermTablePtr dsttable = dst->addTable(name,
                                                      srctable->atomCount(),
                                                      dstparams);
                dsttable->category = srctable->category;

                /* get the terms included in the selection */
                IdList terms = srctable->terms();
                terms.erase(std::remove_if( 
                            terms.begin(), 
                            terms.end(), 
                            BadTerm(srctable, atmmap)),
                        terms.end());
                /* get the set of parameters referenced by at least one term.
                 * The reference counting in ParamTable doesn't help us here
                 * because we aren't interested in params that are referenced
                 * by tables outside of our System.  However, there is a fair
                 * bit of duplication here of logic in append.cxx and that 
                 * should be addressed at some point. */
                IdList params;
                for (Id i=0; i<terms.size(); i++) {
                    Id p = srctable->param(terms[i]);
                    if (!bad(p)) params.push_back(p);
                }
                sort_unique(params);

                /* Add to dstparams just the params that haven't already been
                 * added. */
                IdList newsrcparams;
                for (Id i=0; i<params.size(); i++) {
                    Id p = params[i];
                    if (bad(src2dst.at(p))) newsrcparams.push_back(p);
                }
                IdList newdstparams = AppendParams(dstparams, srcparams, 
                                                   newsrcparams);
                /* update the mapping */
                for (Id i=0; i<newdstparams.size(); i++) {
                    src2dst.at(newsrcparams[i]) = newdstparams[i];
                }

                /* append the terms */
                AppendTerms(dsttable, srctable, atmmap, terms, src2dst);
            }
        }
    }

    /* copy all top-level data */
    dst->name = src->name;
    dst->global_cell = src->global_cell;
    dst->nonbonded_info = src->nonbonded_info;
    for (Provenance const& p : src->provenance()) dst->addProvenance(p);

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

