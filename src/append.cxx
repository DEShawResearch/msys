#include "append.hxx"
#include "system.hxx"
#include "override.hxx"

#include <sstream>
#include <stdexcept>
#include <algorithm>

using namespace desres::msys;

IdList desres::msys::AppendParams( ParamTablePtr dst, 
                                   ParamTablePtr src,
                                   IdList const& params ) {

    IdList ids;

    /* add missing props, and construct a mapping from prop index in 
     * src to prop index in dst */
    Id nprops = src->propCount();
    IdList map(nprops);
    for (Id i=0; i<nprops; i++) {
        map[i] = dst->addProp(src->propName(i), src->propType(i));
    }
    
    /* add params */
    for (Id i=0; i<params.size(); i++) {
        Id id = dst->addParam();
        ids.push_back(id);
        for (Id j=0; j<nprops; j++) {
            dst->value(id,map[j]) = src->value(params[i],j);
        }
    }
    return ids;
}

IdList desres::msys::AppendTerms( TermTablePtr dst, TermTablePtr src, 
                                  IdList const& src2dst,
                                  IdList const& terms ) {

    /* First pass: get the set of parameters that are referenced by at least
     * one term */
    IdList srcparams;
    for (Id i=0; i<terms.size(); i++) {
        Id p = src->param(terms[i]);
        if (!bad(p)) srcparams.push_back(p);
    }
    sort_unique(srcparams);

    /* append the necessary parameters */
    IdList dstparams = AppendParams( dst->params(), src->params(), srcparams );

    /* construct a mapping from srcparams to dstparams */
    IdList pmap(src->params()->paramCount(), BadId);
    for (Id i=0; i<srcparams.size(); i++) pmap[srcparams[i]] = dstparams[i];

    return AppendTerms( dst, src, src2dst, terms, pmap );
}

IdList desres::msys::AppendTerms( TermTablePtr dst, TermTablePtr src, 
                                  IdList const& src2dst,
                                  IdList const& terms,
                                  IdList const& pmap ) {

    IdList ids;

    /* copy over term properties */
    Id nprops = src->termPropCount();
    IdList map(nprops);
    for (Id i=0; i<nprops; i++) {
        map[i] = dst->addTermProp(src->termPropName(i), src->termPropType(i));
    }

    /* append all terms */
    for (Id i=0; i<terms.size(); i++) {
        Id srcterm = terms[i];
        Id srcparam = src->param(srcterm);
        Id dstparam = bad(srcparam) ? BadId : pmap[srcparam];
        IdList atoms = src->atoms(srcterm);
        for (Id j=0; j<atoms.size(); j++) atoms[j] = src2dst[atoms[j]];
        Id dstterm = dst->addTerm(atoms, dstparam);

        for (Id j=0; j<nprops; j++) {
            dst->termPropValue(dstterm,map[j]) = src->termPropValue(srcterm, j);
        }
        ids.push_back(dstterm);
    }

    /* copy overrides */
    if (src->overrides()->count()) {
        IdList dstparams = AppendParams( dst->overrides()->params(),
                                         src->overrides()->params(),
                                         src->overrides()->params()->params());
        std::vector<IdPair> L = src->overrides()->list();
        for (unsigned i=0; i<L.size(); i++) {
            Id p1 = pmap.at(L[i].first);
            Id p2 = pmap.at(L[i].second);
            if (bad(p1) || bad(p2)) continue;
            Id dstparam = dstparams.at(src->overrides()->get(L[i]));
            dst->overrides()->set(IdPair(p1,p2), dstparam);
        }
    }

    return ids;
}

IdList desres::msys::AppendSystem( SystemPtr dstptr, SystemPtr srcptr, Id ctid) {

    System& dst = *dstptr;
    System const& src = *srcptr;

    /* Mappings from src ids to dst ids */
    IdList atmmap(src.maxAtomId(), BadId);
    IdList resmap(src.maxResidueId(), BadId);
    IdList chnmap(src.maxChainId(), BadId);
    IdList ctmap(src.maxCtId(), BadId);

    dst.nonbonded_info.merge(src.nonbonded_info);
    /* only overwrite global cell when appending to the system itself,
     * not to a specific ct. */
    if (bad(ctid)) dst.global_cell.merge(src.global_cell);

    /* copy atom properties */
    Id nprops = src.atomPropCount();
    IdList propmap(nprops);
    for (Id i=0; i<nprops; i++) {
        propmap[i] = dst.addAtomProp(src.atomPropName(i), src.atomPropType(i));
    }

    /* copy bond properties */
    Id nbprops = src.bondPropCount();
    IdList bpropmap(nbprops);
    for (Id i=0; i<nbprops; i++) {
        bpropmap[i] = dst.addBondProp(src.bondPropName(i), src.bondPropType(i));
    }

    /* add cts */
    for (Id srcct : src.cts()) {
        Id dstct = ctid;
        if (bad(dstct)) {
            dstct = dst.addCt();
            dst.ct(dstct) = src.ct(srcct);
        }
        ctmap[srcct] = dstct;
    }

    /* add chains */
    for (Id srcchn : src.chains()) {
        Id srcct = src.chain(srcchn).ct;
        Id dstct = ctmap[srcct];
        Id dstchn = chnmap[srcchn] = dst.addChain(dstct);
        /* copy attributes from src chain to dst chain */
        dst.chain(dstchn) = src.chain(srcchn);
        dst.chain(dstchn).ct = dstct;
    }

    /* add residues */
    for (Id srcres : src.residues()) {
        Id srcchn = src.residue(srcres).chain;
        Id dstchn = chnmap[srcchn];
        Id dstres = resmap[srcres] = dst.addResidue(dstchn);
        /* copy attributes from src residue to dst residue */
        dst.residue(dstres) = src.residue(srcres);
        dst.residue(dstres).chain = dstchn;
    }

    /* add atoms */
    for (Id srcatm : src.atoms()) {
        Id srcres = src.atom(srcatm).residue;
        Id dstres = resmap[srcres];
        Id dstatm = atmmap[srcatm] = dst.addAtom(dstres);
        /* copy attributes from src atom to dst atom */
        dst.atom(dstatm) = src.atom(srcatm);
        dst.atom(dstatm).residue = dstres;
        /* Copy additional atom properties */
        for (Id p=0; p<nprops; p++) {
            dst.atomPropValue(dstatm,propmap[p]) = 
            srcptr->atomPropValue(srcatm, p);
        }
    }

    /* add bonds */
    IdList bonds = src.bonds();
    for (Id i=0; i<bonds.size(); i++) {
        Id srcbnd = bonds[i];
        Id srci = src.bond(srcbnd).i;
        Id srcj = src.bond(srcbnd).j;
        Id dsti = atmmap[srci];
        Id dstj = atmmap[srcj];
        Id dstbnd = dst.addBond(dsti, dstj);
        dst.bond(dstbnd).order = src.bond(srcbnd).order;
        dst.bond(dstbnd).stereo = src.bond(srcbnd).stereo;
        dst.bond(dstbnd).aromatic = src.bond(srcbnd).aromatic;
        dst.bond(dstbnd).resonant_order = src.bond(srcbnd).resonant_order;

        /* Copy additional bond properties */
        for (Id k=0; k<nbprops; k++) {
            dst.bondPropValue(dstbnd,bpropmap[k]) = 
            srcptr->bondPropValue(srcbnd, k);
        }
    }

    /* add/merge term tables */
    std::vector<std::string> tablenames = src.tableNames();
    for (unsigned i=0; i<tablenames.size(); i++) {
        std::string const& name = tablenames[i];
        TermTablePtr srctable = src.table(name);
        TermTablePtr dsttable = dst.table(name);
        if (!dsttable) {
            dsttable = dst.addTable(name, srctable->atomCount());
            dsttable->category = srctable->category;
        } else {
            if (dsttable->category != srctable->category) {
                std::stringstream ss;
                ss << "Append failed: Tables '" << name << "' have different "
                   << "categories: '" << dsttable->category 
                   << "' and '" << srctable->category << "'";
                throw std::runtime_error(ss.str());
            }
        }
        AppendTerms( dsttable, srctable, atmmap, srctable->terms() );
    }

    /* add/replace extra tables */
    std::vector<String> extras = src.auxTableNames();
    for (unsigned i=0; i<extras.size(); i++) {
        std::string const& name = extras[i];
        if (!dst.auxTable(name)) {
            ParamTablePtr srcparams = src.auxTable(name);
            ParamTablePtr dstparams = ParamTable::create();
            AppendParams( dstparams, srcparams, srcparams->params() );
            dst.addAuxTable( name, dstparams );
        }
    }

    return atmmap;
}
