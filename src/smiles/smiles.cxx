#include "smiles.hxx"
#include "../elements.hxx"
#include "../smiles.hxx"

int msys_smiles_parse(desres::msys::smiles::Smiles*);

extern int msys_smiles_debug;

namespace desres { namespace msys { namespace smiles {

    Smiles::Smiles() 
    : txt(), pos(), scanner() 
    {}

    void Smiles::addh(Id atm) {
        Id idH = atoms.size();
        atoms.resize(idH+1);
        atoms[idH].atomic_number = 1;
        Id bnd = bonds.size();
        bonds.resize(bnd+1);
        bonds[bnd].i = atm;
        bonds[bnd].j = idH;
        bonds[bnd].order = 1;
    }

    int Smiles::addh(Id atm, int v1, int v2, int v3) {
        float v = 0;
        for (unsigned i=0; i<bonds.size(); i++) {
            if (bonds[i].i!=atm && bonds[i].j!=atm) continue;
            v += bonds[i].aromatic ? 1.5 : bonds[i].order;
        }
        int h;
        if (v<=v1) {
            h = v1 - v;
        } else if (v<=v2) {
            h = v2 - v;
        } else if (v<=v3) {
            h = v3 - v;
        } else {
            h = 0;
        }
        for (int i=0; i<h; i++) addh(atm);
        return h;
    }

    MoleculePtr Smiles::parse(std::string const& s) {
        if (getenv("MSYS_SMILES_DEBUG")) {
            msys_smiles_debug = 1;
        }

        atoms.clear();
        bonds.clear();
        txt = s.data();
        pos = 0;
        init_scanner();
        int rc = msys_smiles_parse(this);
        msys_smiles_debug = 0;
        destroy_scanner();
        if (rc) {
            std::stringstream ss;
            ss << "parse failed around position " << pos <<":\n" << txt << "\n";
            for (int i=0; i<pos-1; i++) ss << '-';
            ss << "^-\n" << error;
            MSYS_FAIL(ss.str());
        }
        MoleculePtr ptr(new Molecule(atoms.size(), bonds.size()));
        for (unsigned i=0; i<atoms.size(); i++) ptr->atom(i) = atoms[i];
        for (unsigned i=0; i<bonds.size(); i++) ptr->bond(i) = bonds[i];
        return ptr;
    }

    void Smiles::addAtom(atom_t* a, bool organic) {
        char name[4];
        strcpy(name, a->name);
        name[0] = toupper(name[0]);

        a->id = atoms.size();
        atoms.resize(a->id+1);
        atoms[a->id].formal_charge = a->charge;
        atoms[a->id].atomic_number = ElementForAbbreviation(name);
        atoms[a->id].hcount = organic ? -1 : a->hcount;
        atoms[a->id].stereo_parity = a->chiral;
        for (int i=0; i<a->hcount; i++) addh(a->id);
    }

    void Smiles::add(atom_t* a, branch_t* branches) {
        for (branch_t* b = branches; b; b = b->next) {
            add(a, b->chain->first, b->bond);
        }
    }
    void Smiles::add(atom_t* a, ringbond_t* ringbonds) {
        for (ringbond_t* r = ringbonds; r; r = r->next) {
            ringmap::iterator it = rnums.find(r->id);
            if (it == rnums.end()) {
                rnums[r->id] = std::make_pair(a, r->bond);
            } else if (it->second.first->id == a->id) {
                fprintf(stderr, "ringbond %d bonded to itself!", r->id);
            } else if (r->bond != it->second.second &&
                       r->bond != 0 && 
                       it->second.second != 0) {
                fprintf(stderr, "ringbond %d has conflicting bond spec: %c and %c\n", r->id, r->bond, it->second.second);
            } else {
                char bond = r->bond ? r->bond : 
                            it->second.second ? it->second.second :
                            '-';
                add(a, it->second.first, bond);
                rnums.erase(it);
            }
        }
    }

    void Smiles::add(atom_t* ai, atom_t* aj, char bond) {
        if (bond=='.') return;  /* dot means no bond */
        Id id = bonds.size();
        bonds.resize(id+1);
        auto& bnd = bonds[id];
        bnd.i = ai->id;
        bnd.j = aj->id;
        switch (bond) {
            default:
            MSYS_FAIL("Unsupported bond type '" << bond << "'");
            case '-': bnd.order = 1; break;
            case '=': bnd.order = 2; break;
            case '#': bnd.order = 3; break;
        }
        if (bnd.order==1 && islower(ai->name[0]) && islower(aj->name[0])) {
            bnd.aromatic = true;
        }
    }

    void Smiles::finish(chain_t* chain) {
        if (!rnums.empty()) {
            std::stringstream ss;
            ss << "Unclosed rings:\n";
            for (ringmap::iterator it=rnums.begin(); it!=rnums.end(); ++it) {
                ss << it->first << " " << it->second.second << "\n";
            }
            MSYS_FAIL(ss.str());
        }
        /* add hydrogens.  We use the OpenSmiles specification rather than
         * our own AddHydrogen routine, since there could be some differences
         */
        for (unsigned i=0; i<atoms.size(); i++) {
            if (atoms[i].hcount != -1) continue;
            int nh = 0;
            switch (atoms[i].atomic_number) {
                case  5: nh = addh(i, 3); break;
                case  6: nh = addh(i, 4); break;
                case  7: nh = addh(i, 3,5); break;
                case  8: nh = addh(i, 2); break;
                case 15: nh = addh(i, 3,5); break;
                case 16: nh = addh(i, 2,4,6); break;
                case  9:
                case 17:
                case 35:
                case 53: nh = addh(i, 1); break;
                default:;
            }
            atoms[i].hcount = nh;
        }
    }

}}}

namespace desres { namespace msys {

    MoleculePtr FromSmilesString(std::string const& smiles) {
        return smiles::Smiles().parse(smiles);
    }

}}
