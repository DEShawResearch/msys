#include "sitemap.hxx"
#include <stdexcept>
#include <sstream>
#include <cstdio>

using namespace desres::msys::mae;
using desres::fastjson::Json;

SiteMap::SiteMap( SystemPtr h, const Json& sites, const IdList& atoms, 
                  int natoms, int npseudos ) 
: _atoms(atoms) {

    _nsites = sites.get("__size__").as_int();
    if (!_nsites) return;

    const Json& mass = sites.get("ffio_mass");
    const Json& charge = sites.get("ffio_charge");
    const Json& chargeB = sites.get("ffio_chargeB");
    const Json& atype = sites.get("ffio_type");
    const Json& moiety = sites.get("ffio_moiety");
    const Json& vdw = sites.get("ffio_vdwtype");
    const Json& vdwB = sites.get("ffio_vdwtypeB");

    Id iatom = 0;
    Id ipseudo=natoms;
    Id nparticles = natoms + npseudos;
    Id nblocks = nparticles / _nsites;
    _s2p.resize(nparticles);
    for (Id j=0; j<nblocks; j++) {
        for (Id i=0; i<_nsites; i++) {
            Id val;
            if (!strcasecmp(atype.elem(i).as_string("atom"), "atom")) {
                val=iatom++;
            } else {
                val=ipseudo++;
            }
            Id ind = i+j*_nsites;
            _s2p[ind] = val;

            Id atm = atoms[val];
            atom_t& atom = h->atom(atm);

            atom.mass = mass.elem(i).as_float(0);
            const char * vtype = vdw.elem(i).as_string("");
            const char * vtypeB = vdwB.elem(i).as_string(vtype);
            atom.charge  = charge.elem(i).as_float(0);
            atom.chargeB = chargeB.elem(i).as_float(atom.charge);
            atom.moiety  = moiety.elem(i).as_int(0);
            if (chargeB.elem(i).kind()==Json::Float ||
                strcmp(vtype, vtypeB) ||
                atom.moiety!=0) {
                atom.alchemical = true;
            }
        }
    }
}

void SiteMap::addUnrolledTerms(TermTablePtr table, Id param, 
                               const IdList& sites,
                               bool constrained, Id paramB ) const {
    /* do checks on input site ids */
    Id j,m = sites.size();
    for (j=0; j<m; j++) if (sites[j]<1 || sites[j]>_nsites) {
        std::stringstream ss;
        ss << "id " << sites[j] << " is out of bounds (nsites=" << _nsites 
            << ")";
        throw std::runtime_error(ss.str());
    }

    /* find column for marking constrained */
    Id constrained_col = BadId;
    if (constrained) constrained_col = table->termPropIndex("constrained");

    /* unroll the block */
    int i, nblocks = _atoms.size() / _nsites;
    IdList ids(m);
    for (i=0; i<nblocks; i++) {
        int offset = i*_nsites-1;
        /* convert site ids to atom ids */
        for (j=0; j<m; j++) {
            ids[j] = _atoms.at(_s2p.at(sites[j]+offset));
        }
        Id term = table->addTerm(ids, param);
        table->setParamB(term, paramB);
        if (!bad(constrained_col)) {
            table->termPropValue(term, constrained_col)=1;
        }
    }
}

#if 0
Id SiteMap::addTerm( TermTablePtr nb, Id param, const IdList& ids ) const {
    ent::TermTablePtrImpl * nbi = nb.ptr();
    const int m = nb.getAtomsPerItem();
    if (m>(int)ids.size()) FFIO_ERROR("not enough ids supplied");
    for (int j=0; j<m; j++) if (ids[j] < 1 || ids[j] > _nsites ) {
        FFIO_ERROR("id " << ids[j] << " is out of bounds (nsites=" << _nsites);
    }
    std::vector<ent::AtomImpl *> sites(m);
    ent::DictId * de = param.ptr();
    int offset = -1;
    for (int j=0; j<m; j++) {
        sites[j] = _atoms.at(_s2p.at(ids[j]+offset));
    }
    return *nbi->addItem( de, sites );
}
#endif


void SiteMap::addUnrolledPseudoBonds( SystemPtr h, Id parent, Id pseudo) const {
    int i, nblocks = _atoms.size() / _nsites;
    for (i=0; i<nblocks; i++) {
        int offset = i*_nsites-1;
        Id parent_atom = _atoms.at(_s2p.at(parent+offset));
        Id pseudo_atom = _atoms.at(_s2p.at(pseudo+offset));
        h->addBond(parent_atom, pseudo_atom);
        /* asssign pseudo atom to the residue of the parent atom */
        h->setResidue(pseudo_atom, h->atom(parent_atom).residue);
    }
}
