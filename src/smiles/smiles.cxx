#include "smiles.hxx"
#include "../elements.hxx"
#include "../smiles.hxx"
#include <boost/foreach.hpp>

int msys_smiles_parse(desres::msys::smiles::Smiles*);

extern int msys_smiles_debug;

#define ORGANIC AtomProBack
#define INORGANIC AtomOther

namespace desres { namespace msys {
    static void addh(SystemPtr mol, Id atm, int v1, int v2=0, int v3=0) {
        float bonds=0;
        BOOST_FOREACH(Id bnd, mol->bondsForAtom(atm)) {
            bonds += mol->bondFAST(bnd).resonant_order;
        }
        int h;
        if (bonds<=v1) {
            h = v1 - bonds;
        } else if (bonds<=v2) {
            h = v2 - bonds;
        } else if (bonds<=v3) {
            h = v3 - bonds;
        } else {
            h = 0;
        }
        Id res = mol->atomFAST(atm).residue;
        for (int i=0; i<h; i++) {
            Id idH = mol->addAtom(res);
            atom_t& H = mol->atomFAST(idH);
            H.atomic_number = 1;
            H.name = "H";
            mol->addBond(atm,idH);
        }
    }

}}

namespace desres { namespace msys { namespace smiles {

    Smiles::Smiles(SystemPtr m) 
    : txt(), pos(), scanner(), mol(m) {
        Id ct = mol->addCt();
        Id chn = mol->addChain(ct);
        res = mol->addResidue(chn);
    }

    void Smiles::parse(std::string const& s) {
        if (getenv("MSYS_SMILES_DEBUG")) {
            msys_smiles_debug = 1;
        }

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
    }

    void Smiles::addAtom(atom_t* a, bool organic) {
        a->id = mol->addAtom(res);
        msys::atom_t& atm = mol->atomFAST(a->id);
        atm.name = a->name;
        a->name[0] = toupper(a->name[0]);
        atm.atomic_number = ElementForAbbreviation(a->name);
        atm.formal_charge = a->charge;
        atm.type = organic ? ORGANIC : INORGANIC;
        for (int i=0; i<a->hcount; i++) {
            Id idH = mol->addAtom(atm.residue);
            mol->atomFAST(idH).name = "H";
            mol->atomFAST(idH).atomic_number = 1;
            mol->addBond(a->id, idH);
        }
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
                rnums[r->id] = std::make_pair(a->id, r->bond);
            } else if (it->second.first == a->id) {
                fprintf(stderr, "ringbond %d bonded to itself!", r->id);
            } else if (r->bond != it->second.second &&
                       r->bond != 0 && 
                       it->second.second != 0) {
                fprintf(stderr, "ringbond %d has conflicting bond spec: %c and %c\n", r->id, r->bond, it->second.second);
            } else {
                char bond = r->bond ? r->bond : 
                            it->second.second ? it->second.second :
                            '-';
                add(a->id, it->second.first, bond);
                rnums.erase(it);
            }
        }
    }

    void Smiles::add(atom_t* ai, atom_t* aj, char bond) {
        add(ai->id, aj->id, bond);
    }

    void Smiles::add(Id i, Id j, char bond) {
        if (bond=='.') return;  /* dot means no bond */
        Id id = mol->addBond(i,j);
        bond_t& bnd = mol->bondFAST(id);
        switch (bond) {
            default:
            case '-': bnd.order = 1; break;
            case '=': bnd.order = 2; break;
            case '#': bnd.order = 3; break;
            case '$': bnd.order = 4; break;
        }
        if (bnd.order==1 && islower(mol->atomFAST(i).name[0]) &&
                            islower(mol->atomFAST(j).name[0])) {
            bnd.resonant_order = 1.5;
        } else {
            bnd.resonant_order = bnd.order;
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
        for (System::iterator i=mol->atomBegin(), e=mol->atomEnd(); i!=e; ++i) {
            msys::atom_t& atm = mol->atomFAST(*i);
            if (atm.type != ORGANIC) continue;
            int anum = atm.atomic_number;
            switch (anum) {
                case  5: addh(mol, *i, 3); break;
                case  6: addh(mol, *i, 4); break;
                case  7: addh(mol, *i, 3,5); break;
                case  8: addh(mol, *i, 2); break;
                case 15: addh(mol, *i, 3,5); break;
                case 16: addh(mol, *i, 2,4,6); break;
                case  9:
                case 17:
                case 35:
                case 53: addh(mol, *i, 1); break;
                default:;
            }
            mol->atomFAST(*i).type = INORGANIC;
        }
        mol->updateFragids();
    }

}}}

namespace desres { namespace msys {

    SystemPtr FromSmilesString(std::string const& smiles) {
        SystemPtr mol = System::create();
        smiles::Smiles context(mol);
        context.parse(smiles);
        return mol;
    }

}}
