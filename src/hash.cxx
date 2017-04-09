// yes, I'm that evil.  I want to get at the ParamTablePtr objects
//#define private public

#include "hash.hxx"
#include <ThreeRoe/ThreeRoe.hpp>

namespace desres { namespace msys {

    static void hash_params(ThreeRoe& tr, ParamTablePtr params) {
        for (Id iprop=0, nprop=params->propCount(); iprop<nprop; ++iprop) {
            auto name = params->propName(iprop);
            auto type = params->propType(iprop);
            tr.Update(name.data(), name.size());
            tr.Update(&type, sizeof(type));
            for (Id i=0, n=params->paramCount(); i<n; i++) {
                switch (type) {
                case IntType:
                {
                    auto v = params->value(i,iprop).asInt();
                    tr.Update(&v, sizeof(v));
                    break;
                }
                case FloatType:
                {
                    auto v = params->value(i,iprop).asFloat();
                    tr.Update(&v, sizeof(v));
                    break;
                }
                default:
                case StringType:
                {
                    auto v = params->value(i,iprop).c_str();
                    tr.Update(v, strlen(v));
                    break;
                }
                }
            }
        }
    }

    static void hash_atoms(ThreeRoe& tr, SystemPtr mol) {
        for (auto i=mol->atomBegin(), e=mol->atomEnd(); i!=e; ++i) {
            tr.Update(&mol->atomFAST(*i), sizeof(atom_t));
        }
        hash_params(tr, mol->atomProps());
    }

    static void hash_bonds(ThreeRoe& tr, SystemPtr mol) {
        for (auto i=mol->bondBegin(), e=mol->bondEnd(); i!=e; ++i) {
            tr.Update(&mol->bondFAST(*i), sizeof(bond_t));
        }
        hash_params(tr, mol->bondProps());
    }

    static void hash_residues(ThreeRoe& tr, SystemPtr mol) {
        for (auto i=mol->residueBegin(), e=mol->residueEnd(); i!=e; ++i) {
            tr.Update(&mol->residueFAST(*i), sizeof(residue_t));
        }
    }

    static void hash_chains(ThreeRoe& tr, SystemPtr mol) {
        for (auto i=mol->chainBegin(), e=mol->chainEnd(); i!=e; ++i) {
            auto const& chn = mol->chainFAST(*i);
            tr.Update(&chn.ct, sizeof(Id));
            tr.Update(chn.name.data(), chn.name.size());
            tr.Update(chn.segid.data(), chn.segid.size());
        }
    }

    static void hash_cts(ThreeRoe& tr, SystemPtr mol) {
        for (auto ctid : mol->cts()) {
            auto& ct = mol->ctFAST(ctid);
            hash_params(tr, ct.kv());
        }
    }

    namespace {
        struct prop_visitor {
            ThreeRoe& tr;

            prop_visitor(ThreeRoe& t) : tr(t) {}
            void operator()(Int& v) const { tr.Update(&v, sizeof(v)); }
            void operator()(Float& v) const { tr.Update(&v, sizeof(v)); }
            void operator()(String& v) const { tr.Update(v.data(), v.size()); }
        };
    }

    static void hash_tables(ThreeRoe& tr, SystemPtr mol) {
        for (auto name : mol->tableNames()) {
            tr.Update(name.data(), name.size());
            auto table = mol->table(name);
            hash_params(tr, table->params());
            hash_params(tr, table->termProps());
            tr.Update(&table->category, sizeof(table->category));
            for (auto term : *table) {
                tr.Update(term.atoms(), (term.size()+1)*sizeof(Id));
            }
            for (auto& prop : table->tableProps()) {
                tr.Update(prop.first.data(), prop.first.size());
                boost::apply_visitor(prop_visitor(tr), prop.second);
            }

            if (table->overrides()) {
                auto ov = table->overrides();
                hash_params(tr, ov->params());
                for (auto pair : ov->list()) {
                    auto id = ov->get(pair);
                    tr.Update(&id, sizeof(id));
                }

            }
        }
    }

    static void hash_other(ThreeRoe& tr, SystemPtr mol) {
        tr.Update(mol->name.data(), mol->name.size());
        tr.Update(mol->global_cell[0], 9*sizeof(double));
        tr.Update(mol->nonbonded_info.vdw_funct.data(),
                  mol->nonbonded_info.vdw_funct.size());
        tr.Update(mol->nonbonded_info.vdw_rule.data(),
                  mol->nonbonded_info.vdw_rule.size());
        tr.Update(mol->nonbonded_info.es_funct.data(),
                  mol->nonbonded_info.es_funct.size());
        for (auto name : mol->auxTableNames()) {
            tr.Update(name.data(), name.size());
            hash_params(tr, mol->auxTable(name));
        }
    }

    uint64_t HashSystem(SystemPtr mol) {
        ThreeRoe tr;

        hash_atoms(tr, mol);
        hash_bonds(tr, mol);
        hash_residues(tr, mol);
        hash_chains(tr, mol);
        hash_cts(tr, mol);
        hash_tables(tr, mol);
        hash_other(tr, mol);

        return tr.Final().first;
    }
}}
