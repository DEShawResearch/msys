#include "../ff.hxx"
#include "../schema.hxx"
#include <unordered_set>
#include <cmath>

using namespace desres::msys;

struct vdw {
    double sigma;
    double epsilon;
    vdw(double s, double e) : sigma(s), epsilon(e) {}
};

struct lj_pair {
    double aij;
    double bij;
    lj_pair(double a, double b) : aij(a), bij(b) {}
};

typedef vdw (*combining_rule_t)(vdw vi, vdw vj);

static vdw arithmetic_geometric(vdw vi, vdw vj) {
    return vdw(
            0.5 * (vi.sigma + vj.sigma),
            sqrt(vi.epsilon * vj.epsilon));
}

static vdw geometric(vdw vi, vdw vj) {
    return vdw(
            sqrt(vi.sigma * vj.sigma),
            sqrt(vi.epsilon * vj.epsilon));
}

static lj_pair convert_12_6(vdw v, double scale) {
    return lj_pair(
            scale * pow(v.sigma, 12) * v.epsilon * 4.0,
            scale * pow(v.sigma,  6) * v.epsilon * 4.0);
}

static combining_rule_t combiner_for(std::string rule) {
    if (rule == "arithmetic/geometric") return arithmetic_geometric;
    if (rule == "geometric") return geometric;
    MSYS_FAIL("Unsupported combining rule '" << rule << "'");
}

namespace desres { namespace msys { namespace ff {
template<>
void build<Component::exclusions>(
        SystemPtr mol,
        Forcefield const& ff,
        Tuples const& tuples) {

    auto excls = AddTable(mol, "exclusion");
    TermTablePtr pairs_table;
    if (mol->nonbonded_info.vdw_funct == "vdw_12_6") {
        pairs_table = AddTable(mol, "pair_12_6_es");
    } else {
        MSYS_FAIL("Unsupported vdw_funct for pairs: "
                << mol->nonbonded_info.vdw_funct);
    }
    auto nbterms = mol->table("nonbonded");
    auto nbparams = nbterms->params();
    // assuming vdw_12_6 here
    Id sig_propid = nbparams->propIndex("sigma");
    Id eps_propid = nbparams->propIndex("epsilon");
    auto combiner = combiner_for(mol->nonbonded_info.vdw_rule);

    MultiIdList const* all_tuples[3] = {
        &tuples.bonds, &tuples.angles, &tuples.dihedrals
    };
    Id rule = ff.rules.exclusions;
    IdList term(2);

    std::unordered_set<uint64_t> excl_set;

    for (auto tlist : all_tuples) {
        if (tlist->empty()) continue;
        Id arity = tlist->at(0).size();
        if (rule > 0 && arity > rule) continue;
        double es_scale =
            ff.rules.es_scale.empty() ? 0 :
            ff.rules.es_scale.at(arity-2);
        double lj_scale =
            ff.rules.lj_scale.empty() ? 0 :
            ff.rules.lj_scale.at(arity-2);
        for (auto& tuple : *tlist) {
            term[0] = tuple[0];
            term[1] = tuple[arity-1];
            uint64_t key = (uint64_t(term[0]) << 32) | term[1];
            if (!excl_set.emplace(key).second) continue;
            excls->addTerm(term, BadId);
            if (es_scale > 0 || lj_scale > 0) {
                Id ni = nbterms->param(term[0]);
                Id nj = nbterms->param(term[1]);
                double si = nbparams->value(ni, sig_propid);
                double sj = nbparams->value(nj, sig_propid);
                double ei = nbparams->value(ni, eps_propid);
                double ej = nbparams->value(nj, eps_propid);
                auto sij = combiner(vdw(si,ei), vdw(sj,ej));
                auto lj = convert_12_6(sij, lj_scale);

                auto qi = mol->atomFAST(term[0]).charge;
                auto qj = mol->atomFAST(term[1]).charge;
                auto qij = es_scale * (qi * qj);

                auto param = pairs_table->params()->addParam();
                pairs_table->params()->value(param, "aij") = lj.aij;
                pairs_table->params()->value(param, "bij") = lj.bij;
                pairs_table->params()->value(param, "qij") = qij;
                pairs_table->addTerm(term, param);
            }
        }
    }
}

}}} // namespace
