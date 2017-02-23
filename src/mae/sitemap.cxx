#include "sitemap.hxx"
#include <stdexcept>
#include <sstream>
#include <cstdio>

using namespace desres::msys::mae;
using desres::msys::fastjson::Json;

#ifdef WIN32
#define strcasecmp _stricmp
#endif

desres::msys::mae::SiteMap::SiteMap( desres::msys::SystemPtr h, const Json& sites, const IdList& atoms,
                  int natoms, int npseudos ) 
: _atoms(atoms) {

    _nsites = sites.get("__size__").as_int();
    if (!_nsites) return;

    const Json& mass = sites.get("ffio_mass");
    const Json& charge = sites.get("ffio_charge");
    const Json& atype = sites.get("ffio_type");

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
            atom.charge  = charge.elem(i).as_float(0);
        }
    }
}
#ifdef DESMOND_USE_SCHRODINGER_MMSHARE
void desres::msys::mae::SiteMap::addUnrolledTerms(TermTablePtr table, Id param, 
                               const IdList& sites,
                               bool constrained, const char* schedule ) const {
#else
void desres::msys::mae::SiteMap::addUnrolledTerms(TermTablePtr table, Id param, 
                               const IdList& sites,
                               bool constrained ) const {
#endif
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

#ifdef DESMOND_USE_SCHRODINGER_MMSHARE
    /* find column for marking schedule */
    Id schedule_col = BadId;
    if (schedule) schedule_col = table->termPropIndex("schedule");
#endif


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
        if (!bad(constrained_col)) {
            table->termPropValue(term, constrained_col)=1;
        }
#ifdef DESMOND_USE_SCHRODINGER_MMSHARE
        if (!bad(schedule_col)) {
            table->termPropValue(term, schedule_col)=schedule;
        }
#endif
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

#ifdef DESMOND_USE_SCHRODINGER_MMSHARE
void copy_parent_properites(desres::msys::SystemPtr mol,
                            desres::msys::Id parent,
                            desres::msys::Id pseudo)
{

  static const char* ATOM_PROPERTIES[] = { "grp_energy", "grp_ligand" };

  for (unsigned i = 0; i < sizeof(ATOM_PROPERTIES)/sizeof(const char*); i++) { 
    desres::msys::Id pidx = mol->atomPropIndex(ATOM_PROPERTIES[i]);
      if (pidx != desres::msys::BadId) {
        mol->atomPropValue(pseudo, pidx) = mol->atomPropValue(parent, pidx);
      }
  }
}
#endif

void desres::msys::mae::SiteMap::addUnrolledPseudoBonds( desres::msys::SystemPtr h, Id parent, Id pseudo) const {
    int i, nblocks = _atoms.size() / _nsites;
    for (i=0; i<nblocks; i++) {
        int offset = i*_nsites-1;
        Id parent_atom = _atoms.at(_s2p.at(parent+offset));
        Id pseudo_atom = _atoms.at(_s2p.at(pseudo+offset));
        h->addBond(parent_atom, pseudo_atom);

#ifdef DESMOND_USE_SCHRODINGER_MMSHARE
        copy_parent_properites(h, parent_atom, pseudo_atom);
#endif

        /* asssign pseudo atom to the residue of the parent atom */
        h->setResidue(pseudo_atom, h->atom(parent_atom).residue);
    }
}
