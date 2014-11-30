#include "../json.hxx"
#include <fastjson/print.hxx>

using namespace desres::msys;
using desres::msys::fastjson::Json;

namespace {
    void process_particles(SystemPtr sys, Json& js) {
        Json node;
        node.to_object();
        /* built-in properties */
        {
            Json props;
            props.to_array();
            Json prop;
            props.append(prop.to_string("id"));
            props.append(prop.to_string("chain"));
            props.append(prop.to_string("resname"));
            props.append(prop.to_string("resid"));
            props.append(prop.to_string("x"));
            props.append(prop.to_string("y"));
            props.append(prop.to_string("z"));
            props.append(prop.to_string("charge"));
            props.append(prop.to_string("vx"));
            props.append(prop.to_string("vy"));
            props.append(prop.to_string("vz"));
            props.append(prop.to_string("mass"));
            props.append(prop.to_string("chargeB"));
            props.append(prop.to_string("anum"));
            props.append(prop.to_string("formal_charge"));
            node.append("cols", props);
        }
        Json particles;
        particles.to_array();
        for (Id i=0; i<sys->maxAtomId(); i++) {
            if (!sys->hasAtom(i)) continue;
            atom_t const& atm = sys->atom(i);
            residue_t const& res = sys->residue(atm.residue);
            chain_t const& chn = sys->chain(res.chain);
            Json row;
            row.to_array();
            Json val;
            row.append(val.to_int(i));
            row.append(val.to_string(chn.name.c_str()));
            row.append(val.to_string(res.name.c_str()));
            row.append(val.to_int(res.resid));
            row.append(val.to_float(atm.x));
            row.append(val.to_float(atm.y));
            row.append(val.to_float(atm.z));
            row.append(val.to_float(atm.charge));
            row.append(val.to_float(atm.vx));
            row.append(val.to_float(atm.vy));
            row.append(val.to_float(atm.vz));
            row.append(val.to_float(atm.mass));
            row.append(val.to_float(atm.chargeB));
            row.append(val.to_int(atm.atomic_number));
            row.append(val.to_int(atm.formal_charge));
            particles.append(row);
        }
        node.append("rows", particles);
        js.append("particles", node);
    }

    void process_bonds(SystemPtr sys, Json& js) {
        Json node, props, vals, tmp;
        node.to_object();
        props.to_array();
        vals.to_array();

        props.append(tmp.to_string("i"));
        props.append(tmp.to_string("j"));
        props.append(tmp.to_string("order"));

        for (Id i=0; i<sys->maxBondId(); i++) {
            if (!sys->hasBond(i)) continue;
            bond_t const& bnd = sys->bond(i);

            Json row;
            row.to_array();
            row.append(tmp.to_int(bnd.i));
            row.append(tmp.to_int(bnd.j));
            row.append(tmp.to_float(bnd.order));
            vals.append(row);
        }

        node.append("cols", props);
        node.append("rows", vals);
        js.append("bonds", node);
    }
    void process_tables(SystemPtr sys, Json& js) {
        Json node, tmp;
        node.to_object();

        std::vector<String> names = sys->tableNames();
        for (unsigned i=0; i<names.size(); i++) {
            Json tnode;
            tnode.to_object();
            TermTablePtr table = sys->table(names[i]);
            tnode.append("category", tmp.to_string(print(table->category).c_str()));
            tnode.append("natoms", tmp.to_int(table->atomCount()));

            /* term properties */
            Json termprops;
            termprops.to_array();
            for (unsigned j=0; j<table->termPropCount(); j++) {
                termprops.append(tmp.to_string(table->termPropName(j).c_str()));
            }
            tnode.append("termprops", termprops);

            /* param properties */
            ParamTablePtr params = table->params();
            Json props;
            props.to_array();
            for (unsigned j=0; j<params->propCount(); j++) {
                props.append(tmp.to_string(params->propName(j).c_str()));
            }
            tnode.append("props", props);

            Json terms;
            terms.to_array();
            for (Id j=0; j<table->maxTermId(); j++) {
                if (!table->hasTerm(j)) continue;
                Json term;
                term.to_array();
                /* atoms in term */
                Json atoms;
                atoms.to_array();
                for (Id k=0; k<table->atomCount(); k++) {
                    atoms.append(tmp.to_int(table->atom(j,k)));
                }
                term.append(atoms);

                /* term properties */
                termprops.to_array();
                for (unsigned k=0; k<table->termPropCount(); k++) {
                    ValueRef v = table->termPropValue(j,k);
                    switch (v.type()) {
                        case IntType:
                            termprops.append(tmp.to_int(v.asInt()));
                            break;
                        case FloatType:
                            termprops.append(tmp.to_float(v.asFloat()));
                            break;
                        case StringType:
                        default:
                            termprops.append(tmp.to_string(v.asString().c_str()));
                    }
                }
                term.append(termprops);

                /* param properties */
                props.to_array();
                Id param = table->param(j);
                for (unsigned k=0; k<params->propCount(); k++) {
                    ValueRef v = params->value(param,k);
                    switch (v.type()) {
                        case IntType:
                            props.append(tmp.to_int(v.asInt()));
                            break;
                        case FloatType:
                            props.append(tmp.to_float(v.asFloat()));
                            break;
                        case StringType:
                        default:
                            props.append(tmp.to_string(v.asString().c_str()));
                    }
                }
                term.append(props);

                /* alchemical param properties */
                Id paramB = table->paramB(j);
                if (bad(paramB)) {
                    term.append(props.to_null());
                } else {
                }
                terms.append(term);
            }
            tnode.append("terms", terms);
            node.append(names[i].c_str(), tnode);
        }
        js.append("tables", node);
    }
}

void desres::msys::ExportJSON( SystemPtr sys, std::ostream& out ) {
    Json js;
    js.to_object();
    process_particles(sys,js);
    process_bonds(sys,js);
    process_tables(sys,js);
    print_json(out, js);
}
