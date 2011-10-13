#include "../json.hxx"
#include <fastjson/print.hxx>

using namespace desres::msys;
using desres::fastjson::Json;
using desres::fastjson::Json;

namespace {
    void process_particles(SystemPtr sys, Json& js) {
        Json node;
        node.to_object();
        /* built-in properties */
        {
            Json props;
            props.to_array();
            Json prop;
            props.append(prop.to_string("gid"));
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
            row.append(val.to_int(atm.gid));
            row.append(val.to_string(chn.name.c_str()));
            row.append(val.to_string(res.name.c_str()));
            row.append(val.to_int(res.num));
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
            atom_t const& ai = sys->atom(bnd.i);
            atom_t const& aj = sys->atom(bnd.j);

            Json row;
            row.to_array();
            row.append(tmp.to_int(ai.gid));
            row.append(tmp.to_int(aj.gid));
            row.append(tmp.to_float(bnd.order));
            vals.append(row);
        }

        node.append("cols", props);
        node.append("rows", vals);
        js.append("bonds", node);
    }
}

void desres::msys::ExportJSON( SystemPtr sys, std::ostream& out ) {
    Json js;
    js.to_object();
    process_particles(sys,js);
    process_bonds(sys,js);
    print_json(out, js);
}
