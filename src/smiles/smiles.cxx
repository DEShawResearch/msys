#include "smiles.hxx"
#include "../elements.hxx"
#include "../smiles.hxx"

int msys_smiles_parse(desres::msys::smiles::Smiles*);

extern int msys_smiles_debug;

namespace desres { namespace msys { namespace smiles {

    Smiles::Smiles() 
    : txt(), pos(), scanner() 
    {}

    void Smiles::addh(Id id) {
        auto idH = mol->addAtom(0);
        auto& atm = mol->atomFAST(idH);
        atm.atomic_number = 1;
        atm.name = "H";
        mol->addBond(id, idH);
    }

    int Smiles::addh(Id atm, int v1, int v2, int v3) {
        float v = 0;
        for (unsigned i=0; i<mol->maxBondId(); i++) {
            auto const& bnd = mol->bondFAST(i);
            if (bnd.i != atm && bnd.j != atm) continue;
            v += bnd.aromatic ? 1.5 : bnd.order;
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

    SystemPtr Smiles::parse(std::string const& s) {
        if (getenv("MSYS_SMILES_DEBUG")) {
            msys_smiles_debug = 1;
        }

        hcount.clear();
        mol = System::create();
        mol->addChain();
        mol->addResidue(0);
        mol->ct(0).setName(s);
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
        return mol;
    }

    void Smiles::addAtom(atom_t* a, bool organic) {
        char name[4];
        strcpy(name, a->name);
        name[0] = toupper(name[0]);

        a->id = mol->addAtom(0);
        auto& atm = mol->atomFAST(a->id);
        atm.formal_charge = a->charge;
        atm.atomic_number = ElementForAbbreviation(name);
        hcount.push_back(organic ? -1 : a->hcount);
        atm.stereo_parity = a->chiral;
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
        Id id = mol->addBond(ai->id, aj->id);
        auto& bnd = mol->bondFAST(id);
        switch (bond) {
            default:
            MSYS_FAIL("Unsupported bond type '" << bond << "'");
            case '/':
            case '\\':
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
        for (unsigned i=0; i<mol->maxAtomId(); i++) {
            if (hcount[i] != -1) continue;
            switch (mol->atomFAST(i).atomic_number) {
                case  5: addh(i, 3); break;
                case  6: addh(i, 4); break;
                case  7: addh(i, 3,5); break;
                case  8: addh(i, 2); break;
                case 15: addh(i, 3,5); break;
                case 16: addh(i, 2,4,6); break;
                case  9:
                case 17:
                case 35:
                case 53: addh(i, 1); break;
                default:;
            }
        }
    }

}}}

namespace desres { namespace msys {

    SystemPtr FromSmilesString(std::string const& smiles) {
        return smiles::Smiles().parse(smiles);
    }

}}
