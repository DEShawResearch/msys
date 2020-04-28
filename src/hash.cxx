#include "hash.hxx"
#include "MsysThreeRoe.hpp"

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
            // atoms memset 0 themselves, so this is ok.
            tr.Update(&mol->atomFAST(*i), sizeof(atom_t));
        }
        hash_params(tr, mol->atomProps());
    }

    static void hash_bonds(ThreeRoe& tr, SystemPtr mol) {
        for (auto i=mol->bondBegin(), e=mol->bondEnd(); i!=e; ++i) {
            // bonds don't memset zero themselves, so we must be careful.
            auto& b = mol->bondFAST(*i);
            tr.Update(&b.i, sizeof(b.i));
            tr.Update(&b.j, sizeof(b.j));
            tr.Update(&b.order, sizeof(b.order));
            tr.Update(&b.stereo, sizeof(b.stereo));
            tr.Update(&b.aromatic, sizeof(b.aromatic));
        }
        hash_params(tr, mol->bondProps());
    }

    static void hash_residues(ThreeRoe& tr, SystemPtr mol) {
        for (auto i=mol->residueBegin(), e=mol->residueEnd(); i!=e; ++i) {
            // even though this struct is POD, smallstrings don't bzero
            auto& r = mol->residueFAST(*i);
            tr.Update(&r.chain,sizeof(r.chain));
            tr.Update(&r.resid,sizeof(r.resid));
            tr.Update(r.name.c_str(),r.name.size());
            tr.Update(r.insertion.c_str(),r.insertion.size());
            tr.Update(&r.resid,sizeof(r.resid));
            tr.Update(&r.type,sizeof(r.type));
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
            using result_type = void;
            ThreeRoe& tr;

            prop_visitor(ThreeRoe& t) : tr(t) {}
            void operator()(Int const& v) const { tr.Update(&v, sizeof(v)); }
            void operator()(Float const& v) const { tr.Update(&v, sizeof(v)); }
            void operator()(String const& v) const { tr.Update(v.data(), v.size()); }
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
        // don't include system name
        //tr.Update(mol->name.data(), mol->name.size());
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
        static const bool debug=false;

        if (debug) std::cout << tr.Final().first << '\n';
        hash_atoms(tr, mol);
        if (debug) std::cout << tr.Final().first << '\n';
        hash_bonds(tr, mol);
        if (debug) std::cout << tr.Final().first << '\n';
        hash_residues(tr, mol);
        if (debug) std::cout << tr.Final().first << '\n';
        hash_chains(tr, mol);
        if (debug) std::cout << tr.Final().first << '\n';
        hash_cts(tr, mol);
        if (debug) std::cout << tr.Final().first << '\n';
        hash_tables(tr, mol);
        if (debug) std::cout << tr.Final().first << '\n';
        hash_other(tr, mol);
        if (debug) std::cout << tr.Final().first << '\n';

        return tr.Final().first;
    }
}}

