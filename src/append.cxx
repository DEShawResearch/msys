#include "append.hxx"
#include "system.hxx"

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
        Id pB = src->paramB(terms[i]);
        if (!bad(pB)) srcparams.push_back(pB);
    }
    std::sort(srcparams.begin(), srcparams.end());
    srcparams.resize(
            std::unique(srcparams.begin(), srcparams.end())-srcparams.begin());

    /* append the necessary parameters */
    IdList dstparams = AppendParams( dst->params(), src->params(), srcparams );

    /* construct a mapping from srcparams to dstparams */
    IdList pmap(src->params()->paramCount(), BadId);
    for (Id i=0; i<srcparams.size(); i++) pmap[srcparams[i]] = dstparams[i];

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
        Id srcparamB = src->paramB(srcterm);
        Id dstparam = bad(srcparam) ? BadId : pmap[srcparam];
        Id dstparamB = bad(srcparamB) ? BadId : pmap[srcparamB];
        IdList atoms = src->atoms(srcterm);
        for (Id j=0; j<atoms.size(); j++) atoms[j] = src2dst[atoms[j]];
        Id dstterm = dst->addTerm(atoms, dstparam);
        dst->setParamB(dstterm, dstparamB);

        for (Id j=0; j<nprops; j++) {
            dst->termPropValue(dstterm,map[j]) = src->termPropValue(srcterm, j);
        }
        ids.push_back(dstterm);
    }
    return ids;
}

IdList desres::msys::AppendSystem( SystemPtr dstptr, SystemPtr srcptr ) {

    /* Disallow if the the systems have different nonbonded vdw_funct */
    String const& dstvdw = dstptr->nonbonded_info.vdw_funct;
    String const& srcvdw = srcptr->nonbonded_info.vdw_funct;
    if (dstvdw.size() && srcvdw.size() && dstvdw != srcvdw) {
        std::stringstream ss;
        ss << "AppendSystem: incompatible vdw_funct: '" << dstvdw 
           << "' != '" << srcvdw << "'";
        throw std::runtime_error(ss.str());
    }

    System& dst = *dstptr;
    System const& src = *srcptr;
    IdList src2dst(src.maxAtomId(), BadId);

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

    /* add chains */
    IdList chains = src.chains();
    for (Id i=0; i<chains.size(); i++) {
        Id srcchn = chains[i];
        Id dstchn = dst.addChain();
        /* copy attributes from src chain to dst chain */
        dst.chain(dstchn) = src.chain(srcchn);

        /* add residues */
        IdList const& residues = src.residuesForChain(srcchn);
        for (Id j=0; j<residues.size(); j++) {
            Id srcres = residues[j];
            Id dstres = dst.addResidue(dstchn);
            /* copy attributes from src residue to dst residue */
            dst.residue(dstres) = src.residue(srcres);
            dst.residue(dstres).chain = dstchn;

            /* add atoms */
            IdList const& atoms = src.atomsForResidue(srcres);
            for (Id k=0; k<atoms.size(); k++) {
                Id srcatm = atoms[k];
                Id dstatm = dst.addAtom(dstres);
                /* copy attributes from src atom to dst atom */
                dst.atom(dstatm) = src.atom(srcatm);
                dst.atom(dstatm).residue = dstres;
                /* Copy additional atom properties */
                for (Id p=0; p<nprops; p++) {
                    dst.atomPropValue(dstatm,propmap[p]) = 
                    srcptr->atomPropValue(srcatm, p);
                }
                /* map src to dst atoms so we can add bonds */
                src2dst[srcatm] = dstatm;
            }
        }
    }
    /* add bonds */
    IdList bonds = src.bonds();
    for (Id i=0; i<bonds.size(); i++) {
        Id srcbnd = bonds[i];
        Id srci = src.bond(srcbnd).i;
        Id srcj = src.bond(srcbnd).j;
        Id dsti = src2dst[srci];
        Id dstj = src2dst[srcj];
        Id dstbnd = dst.addBond(dsti, dstj);
        dst.bond(dstbnd).order = src.bond(srcbnd).order;

        /* Copy additional bond properties */
        for (Id k=0; k<nbprops; k++) {
            dst.bondPropValue(dstbnd,bpropmap[k]) = 
            srcptr->bondPropValue(srcbnd, k);
        }
    }

    /* add glue */
    std::vector<glue_t> glue = src.gluePairs();
    for (Id i=0; i<glue.size(); i++) {
        Id dsti = src2dst.at(glue[i].first);
        Id dstj = src2dst.at(glue[i].second);
        dst.addGluePair(dsti, dstj);
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
        AppendTerms( dsttable, srctable, src2dst, srctable->terms() );
    }

    /* add/replace extra tables */
    std::vector<String> extras = src.auxTableNames();
    for (unsigned i=0; i<extras.size(); i++) {
        std::string const& name = extras[i];
        ParamTablePtr srcparams = src.auxTable(name);
        ParamTablePtr dstparams = ParamTable::create();
        AppendParams( dstparams, srcparams, srcparams->params() );
        dst.addAuxTable( name, dstparams );
    }

    return src2dst;
}
