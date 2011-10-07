/* @COPYRIGHT@ */

#include "../mae.hxx"
#include "../schema.hxx"
#include "../term_table.hxx"
#include "../clone.hxx"
#include "destro/destro/Destro.hxx"

#include <fstream>
#include <string>
#include <sstream>
#include <stdexcept>
#include <cmath>

using desres::Maeff;
using desres::Destro;
using desres::DestroArray;

using namespace desres::msys;

static std::string pad( const std::string &s, unsigned pad ) {
    std::string result(s);
    static const std::string space(" ");

    for (unsigned i=s.size(); i<pad; i++) {
        if (i%2) result = result + space;
        else     result = space + result;
    }

    return result;
}

static void build_ct_fields( SystemPtr mol, Destro& M ) {
    /* global cell */
    for (int i=0; i<3; i++) {
        for (int j=0; j<3; j++) {
            char schema[32];
            sprintf(schema, "chorus_box_%c%c", 'a'+i, 'x'+j);
            M.add_schema( 'r', schema );
            M[schema] = mol->global_cell[i][j];
        }
    }
}

static void build_m_atom( SystemPtr mol, Destro& M ) {
    /* the fields for the m_atom array come from the VMD maeff plugin */
    static const char * fields[] = {
        "i_m_mmod_type",
        "r_m_x_coord",
        "r_m_y_coord",
        "r_m_z_coord",
        "i_m_residue_number",
        "s_m_insertion_code",
        "s_m_mmod_res",
        "s_m_chain_name",
        "i_m_color",
        "r_m_charge1",
        "r_m_charge2",
        "s_m_pdb_residue_name",
        "s_m_pdb_atom_name", 
        "s_m_grow_name",
        "i_m_atomic_number",
        "i_m_visibility",
        "s_m_pdb_segment_name",
        "r_ffio_x_vel",
        "r_ffio_y_vel",
        "r_ffio_z_vel",
        "i_m_formal_charge"
    };
    static const std::string blank(" ");


    DestroArray& m_atom = M.new_array( "m_atom" );
    for (unsigned i=0; i<sizeof(fields)/sizeof(fields[0]); i++) {
        m_atom.add_schema( fields[i][0], fields[i]+2 );
    }

    IdList ids=mol->atoms();
    for (Id i=0; i<ids.size(); i++) {
        Id id=ids[i];
        const atom_t& atm = mol->atom(id);
        int color, mmod;
        if (atm.atomic_number<1) continue;
        Destro& rec = m_atom.append();
        switch  (atm.atomic_number) {
            default: color=2;  mmod=64; break;  // gray; "any atom"
            case 1:  color=21; mmod=48; break;  // H
            case 3:  color=4;  mmod=11; break;  // Li+ ion
            case 6:  color=2 ; mmod=14; break;  // C
            case 7:  color=43; mmod=40; break;  // N
            case 8:  color=70; mmod=23; break;  // O
            case 9:  color=8;  mmod=56; break;  // F
            case 11: color=4;  mmod=66; break;  // Na+ ion
            case 12: color=4;  mmod=72; break;  // Mg2+ ion
            case 14: color=14; mmod=60; break;  // Si
            case 15: color=15; mmod=53; break;  // P
            case 16: color=13; mmod=52; break;  // S
            case 17: color=13; mmod=102; break;  // Cl- ion
            case 19: color=4;  mmod=67; break;  // K+ ion
            case 20: color=4;  mmod=70; break;  // Ca2+ ion
        }
        const residue_t& res = mol->residue(atm.residue);
        const chain_t& chn = mol->chain(res.chain);

        rec["m_pdb_atom_name"] = pad(atm.name, 4);
        rec["m_pdb_residue_name"] = pad(res.name, 4);
        rec["m_chain_name"] = chn.name;
        rec["m_pdb_segment_name"] = chn.name;
        rec["m_residue_number"] = res.num;
        rec["m_x_coord"] = atm.x;
        rec["m_y_coord"] = atm.y;
        rec["m_z_coord"] = atm.z;
        rec["ffio_x_vel"] = atm.vx;
        rec["ffio_y_vel"] = atm.vy;
        rec["ffio_z_vel"] = atm.vz;
        rec["m_atomic_number"] = atm.atomic_number;
        rec["m_mmod_type"] = mmod;
        rec["m_color"] = color;
        rec["m_visibility"] = 1;
        rec["m_charge1"] = 0.;
        rec["m_charge2"] = 0.;
        rec["m_mmod_res"] = blank;
        rec["m_grow_name"] = blank;
        rec["m_insertion_code"] = " ";
        rec["m_formal_charge"] = atm.formal_charge;

        static const char * entprops[] = {
            "grp_temperature", 
            "grp_energy", 
            "grp_frozen", 
            "grp_bias", 
            "grp_ligand" };
        static const char * maeprops[] = {
            "ffio_grp_thermostat",
            "ffio_grp_energy",
            "ffio_grp_frozen",
            "ffio_grp_cm_moi",
            "ffio_grp_ligand" };
        for (unsigned j=0; j<sizeof(entprops)/sizeof(entprops[0]); j++) {
            Id col = mol->atomPropIndex(entprops[j]);
            if (!bad(col)) {
                rec.add_schema('i', maeprops[j]);
                rec[maeprops[j]]=mol->atomPropValue(id,col).asInt();
            }
        }
    }
}

static void build_m_bond( SystemPtr mol, Destro& M ) {
    DestroArray& m_bond = M.new_array("m_bond");
    m_bond.add_schema('i', "m_from");
    m_bond.add_schema('i', "m_to");
    m_bond.add_schema('i', "m_order");

    IdList ids=mol->bonds();
    for (Id i=0; i<ids.size(); i++) {
        const bond_t& bond = mol->bond(ids[i]);
        Destro& row = m_bond.append();
        row["m_from"] = 1+bond.i;
        row["m_to"]   = 1+bond.j;
        row["m_order"] = (int)bond.order;
    }
}

static void build_sites( SystemPtr mol, Destro& M) {
    static const std::string ATOM("atom");
    static const std::string PSEUDO("pseudo");
    IdList ids=mol->atoms();
    DestroArray& ffio_sites = M.new_array("ffio_sites");
    ffio_sites.add_schema('s', "ffio_type");
    ffio_sites.add_schema('r', "ffio_charge");
    ffio_sites.add_schema('r', "ffio_mass");

    for (unsigned i=0; i<ids.size(); i++) {
        Id id=ids[i];
        const atom_t& atm = mol->atom(id);
        Destro& row = ffio_sites.append();
        row["ffio_type"] = atm.atomic_number==0 ? PSEUDO : ATOM;
        row["ffio_charge"] = atm.charge;
        row["ffio_mass"] = atm.mass;
    }
}

static void build_pseudos( SystemPtr mol, Destro& M ) {
    static const char * schema[] = {
        "r_ffio_x_coord",
        "r_ffio_y_coord",
        "r_ffio_z_coord",
        "s_ffio_atom_name",
        "s_ffio_pdb_residue_name",
        "s_ffio_chain_name",
        "s_ffio_pdb_segment_name",
        "i_ffio_residue_number",
        "r_ffio_x_vel",
        "r_ffio_y_vel",
        "r_ffio_z_vel"
    };
    DestroArray& pseudos = M.new_array("ffio_pseudo");
    for (unsigned i=0; i<sizeof(schema)/sizeof(schema[0]); i++) {
        pseudos.add_schema(schema[i][0], schema[i]+2);
    }

    IdList ids=mol->atoms();
    for (unsigned i=0; i<ids.size(); i++) {
        Id id=ids[i];
        const atom_t& atm = mol->atom(id);
        const residue_t& res = mol->residue(atm.residue);
        const chain_t& chn = mol->chain(res.chain);
        if (atm.atomic_number!=0) continue;
        Destro& row = pseudos.append();
        row["ffio_atom_name"] = atm.name;
        row["ffio_x_coord"] = atm.x;
        row["ffio_y_coord"] = atm.y;
        row["ffio_z_coord"] = atm.z;
        row["ffio_x_vel"] = atm.vx;
        row["ffio_y_vel"] = atm.vy;
        row["ffio_z_vel"] = atm.vz;
        row["ffio_pdb_residue_name"] = pad(res.name,4);
        row["ffio_chain_name"] = chn.name;
        row["ffio_pdb_segment_name"] = chn.name;
        row["ffio_residue_number"] = res.num;
    }
}

struct dms_to_mae {
    const char * dms;     /* table name in dms file */
    const char * mae;     /* ffio_ff subblock name */
    const char * funct;   /* ffio_funct string */
    int nsites;           /* number of automatic p0, p1, .. conversions */
    const char ** params; /* NULL-terminated list mapping to ffio_c1, ... */

    /* This function handles everything not automatically handled by nsites
     * and params.  For example, if a column in the dms has to be mapped to
     * ffio_c0 instead of c1, you'd handle it here.  Similarly, since the
     * virtual sites have a funny ffio_index instead of ffio_ai field,
     * the handler for virtuals lists nsites=0 and handles the site mapping
     * manually.  */
    void (*apply)(TermTablePtr table, Id term, Destro& row);
};


static void constrained_apply(TermTablePtr table, Id term, Destro& row) {
    Id col;
    if (table->findTermProp("constrained", &col) &&
        table->termPropValue(term, col).asInt()) {
        row["ffio_funct"]=std::string(row["ffio_funct"])+"_constrained";
    }
}

static const char * stretch_harm_params[] = { "r0", "fc", NULL };
static const char * angle_harm_params[] = { "theta0", "fc", NULL };
static const char * dihedral_trig_params[] = {
    "fc0", "fc1", "fc2", "fc3", "fc4", "fc5", "fc6", NULL
};
static const char * improper_harm_params[] = { "fc", NULL };

static void dihedral_trig_apply(TermTablePtr table, Id term, Destro& row) {
    row.add_schema('r', "ffio_c0");
    row["ffio_c0"] = table->propValue(term,"phi0").asFloat();
}

static const char * pair_charge_params[] = { "qij", NULL };
static void pair_lj12_6_apply( TermTablePtr table, Id term, Destro& row) {
    row.add_schema('r', "ffio_c1");
    row.add_schema('r', "ffio_c2");
    double aij = table->propValue(term,"aij").asFloat();
    double bij = table->propValue(term,"bij").asFloat();
    double sij=1, eij=0;
    if (aij!=0 && bij!=0) {
        sij = pow(aij/bij, 1./6.);
        eij = (bij*bij) / (4*aij);
    }
    row["ffio_c1"] = sij;
    row["ffio_c2"] = eij;
}

static void cmap_apply( TermTablePtr table, Id term, Destro& row) {
    row.add_schema( 'i', "ffio_c1");
    std::string cmapid = table->propValue(term, "cmapid").asString();
    int id;
    if (sscanf(cmapid.c_str(), "cmap%d", &id)!=1) {
        std::stringstream ss;
        ss << "Invalid cmapid '" << cmapid << "'";
        throw std::runtime_error(ss.str());
    }
    row["ffio_c1"] = id;
}

static void build_extra( SystemPtr mol, Destro& ff ) {
    std::vector<std::string> extras = mol->extraNames();
    for (unsigned i=0; i<extras.size(); i++) {
        const std::string& name = extras[i];
        if (name.substr(0,4)=="cmap") {
            std::string blockname("ffio_cmap");
            blockname += name.substr(4);
            printf("addding cmap table %s\n", blockname.c_str());
            DestroArray& cmap = ff.new_array(blockname);
            cmap.add_schema('r', "ffio_ai");
            cmap.add_schema('r', "ffio_aj");
            cmap.add_schema('r', "ffio_c1");

            ParamTablePtr d = mol->extra(name);
            unsigned phicol = d->propIndex("phi");
            unsigned psicol = d->propIndex("psi");
            unsigned enecol = d->propIndex("energy");
            for (unsigned i=0; i<d->paramCount(); i++) {
                Destro& row = cmap.append();
                row["ffio_ai"] = d->value(i,phicol).asFloat();
                row["ffio_aj"] = d->value(i,psicol).asFloat();
                row["ffio_c1"] = d->value(i,enecol).asFloat();
            }
        }
    }
}

static const char * ah1_params[] = { "r1", NULL };
static const char * ah2_params[] = { "r1", "r2", NULL };
static const char * ah3_params[] = { "r1", "r2", "r3", NULL };
static const char * ah4_params[] = { "r1", "r2", "r3", "r4", NULL };
static const char * hoh_params[] = { "theta", "r1", "r2", NULL };

static const char * lc2_params[] = { "c1", NULL };
static const char * lc3_params[] = { "c1", "c2", NULL };
static const char * out3_params[] = { "c1", "c2", "c3", NULL };
static void virtual_apply( TermTablePtr table, Id term, Destro& row) {
    IdList atoms = table->atoms(term);
    row.add_schema('i', "ffio_index");
    row["ffio_index"] = 1+atoms[0];
    for (unsigned i=1; i<atoms.size(); i++) {
        char buf[32];
        sprintf(buf, "ffio_a%c", 'i'+i-1);
        row.add_schema('i', buf);
        row[buf] = 1+atoms[i];
    }
}

static const char * posre_params[] = { "fcx", "fcy", "fcz", NULL };
static void posre_apply( TermTablePtr table, Id term, Destro& row) {
    row.add_schema('r', "ffio_t1");
    row.add_schema('r', "ffio_t2");
    row.add_schema('r', "ffio_t3");
    /* we've been putting the x0, y0, z0 in the posre_harm_param table,
     * but really they should be item props since they'll be different for
     * every item.  Check both places. */
    const char * props[] = { "x0", "y0", "z0" };
    const char * cols[] = { "ffio_t1", "ffio_t2", "ffio_t3" };
    for (int i=0; i<3; i++) {
        const char * prop = props[i];
        const char * col = cols[i];
        if (table->findTermProp(prop)) {
            row[col]=table->termPropValue(term,prop).asFloat();
        } else {
            row[col]=table->propValue(term,prop).asFloat();
        }
    }
}

static const char * improper_anharm_params[] = { "fc2", "fc4", NULL };
static const char * inplanewag_params[] = { "w0", "fc", NULL };
static const char * pseudopol_params[] = { "a", "b", "cutoff", NULL };


static dms_to_mae dtm_map[] = {
    { "stretch_harm", "ffio_bonds", "harm", 2, stretch_harm_params, constrained_apply },
    { "angle_harm", "ffio_angles", "harm", 3, angle_harm_params, constrained_apply },
    { "dihedral_trig", "ffio_dihedrals", "proper_trig", 
        4, dihedral_trig_params, dihedral_trig_apply },
    { "improper_harm", "ffio_dihedrals", "improper_harm",
        4, improper_harm_params, dihedral_trig_apply },
    { "pair_12_6_es", "ffio_pairs", "coulomb_qij",
        2, pair_charge_params },
    { "pair_12_6_es", "ffio_pairs", "lj12_6_sig_epsilon",
        2, NULL, pair_lj12_6_apply },
    { "torsiontorsion_cmap", "ffio_torsion_torsion", "cmap", 8, 
        NULL, cmap_apply },
    { "constraint_ah1", "ffio_constraints", "ah1", 2, ah1_params },
    { "constraint_ah2", "ffio_constraints", "ah2", 3, ah2_params },
    { "constraint_ah3", "ffio_constraints", "ah3", 4, ah3_params },
    { "constraint_ah4", "ffio_constraints", "ah4", 5, ah4_params },
    { "constraint_hoh", "ffio_constraints", "hoh", 3, hoh_params },
    { "virtual_lc2", "ffio_virtuals", "lc2", 0, lc2_params, virtual_apply },
    { "virtual_lc3", "ffio_virtuals", "lc3", 0, lc3_params, virtual_apply },
    { "virtual_out3", "ffio_virtuals", "out3", 0, out3_params, virtual_apply },
    { "posre_harm", "ffio_restraints", "harm", 1, posre_params, posre_apply },
    { "exclusion", "ffio_exclusions", NULL, 2 },
    { "improper_anharm", "ffio_dihedrals", "improper_anharm", 4, improper_anharm_params },
    { "inplanewag_harm", "ffio_inplanewags", "harm", 4, inplanewag_params },
    { "pseudopol_fermi", "ffio_pseudo_polarization", "fermi", 4, pseudopol_params },
};


static 
void build_tuple_table( SystemPtr mol, TermTablePtr table, 
                        const std::string& name, Destro& ffio_ff ) {

    int handled = false;
    for (unsigned i=0; i<sizeof(dtm_map)/sizeof(dtm_map[0]); i++) {
        const dms_to_mae& dtm = dtm_map[i];
        if (name!=dtm.dms) continue;
        handled = true;

        /* build subblock if we haven't already */
        DestroArray * arr;
        if (ffio_ff.has_block(dtm.mae)) {
            arr = dynamic_cast<DestroArray*>(&ffio_ff.block(dtm.mae));
        } else {
            arr = &ffio_ff.new_array(dtm.mae);
        }
        /* add schema for sites */
        for (int j=0; j<dtm.nsites; j++) {
            char buf[32];
            sprintf(buf, "ffio_a%c", 'i'+j);
            arr->add_schema('i', buf);
        }
        /* add schema for funct */
        if (dtm.funct) arr->add_schema('s', "ffio_funct");
        /* add schema for params */
        if (dtm.params) for (int j=0; dtm.params[j]; j++) {
            char buf[32];
            sprintf(buf, "ffio_c%d", j+1);
            arr->add_schema('r', buf);
        }

        /* iterate over terms */
        IdList terms = table->terms();
        for (unsigned i=0; i<terms.size(); i++) {
            Destro& row = arr->append();
            Id term = terms[i];
            /* add sites */
            IdList atoms = table->atoms(term);
            for (int j=0; j<dtm.nsites; j++) {
                char to[32];
                sprintf(to, "ffio_a%c", 'i'+j);
                row[to] = 1+atoms[j];
            }
            /* add funct */
            if (dtm.funct) row["ffio_funct"] = dtm.funct;
            /* add params */
            if (dtm.params) {
                Id param = table->param(term);
                for (int j=0; dtm.params[j]; j++) {
                    Id col = table->propIndex(dtm.params[j]);
                    char to[32];
                    sprintf(to, "ffio_c%d", j+1);
                    row[to]=table->propValue(param,col).asFloat();
                }
            }
            if (dtm.apply) dtm.apply( table, term, row );
        }
    }
    if (!handled) {
        std::stringstream ss;
        ss << "Failed to process dms table '" << name << "'";
        std::cerr << ss.str() << std::endl;
        //throw std::runtime_error(ss.str());
    }
}

static void build_nonbonded(SystemPtr mol, TermTablePtr table, Destro& ffio_ff) {
    /* grab nonbonded info for combining rule and funct */
    std::string funct = mol->nonbonded_info.vdw_funct;
    std::string rule = mol->nonbonded_info.vdw_rule;

    /* translate the funct from dms to mae */
    static const char * vdw_12_6[] = { "sigma", "epsilon" };
    static const char * vdw_exp_6[] = { "alpha", "epsilon", "rmin" };
    const char ** vdwprops = NULL;
    int nprops=0;
    if (funct=="vdw_12_6") {
        funct="LJ12_6_sig_epsilon";
        vdwprops = vdw_12_6;
        nprops=2;
    }
    else if (funct=="vdw_exp_6") {
        funct="exp_6x";
        vdwprops = vdw_exp_6;
        nprops=3;
    } else {
        std::stringstream ss;
        ss << "Unsupported vdw_funct '" << ss.str() << "'";
        throw std::runtime_error(ss.str());
    }

    /* set the combining rule in the ffio_ff block */
    ffio_ff.add_schema('s', "ffio_comb_rule");
    ffio_ff["ffio_comb_rule"] = rule;

    /* Add the vdwtype schema to ffio_sites */
    Destro &sites = ffio_ff.block("ffio_sites");
    sites.add_schema('s', "ffio_vdwtype");

    /* construct the vdwtypes schema */
    DestroArray& vdw = ffio_ff.new_array("ffio_vdwtypes");
    vdw.add_schema('s', "ffio_name");
    vdw.add_schema('s', "ffio_funct");

    /* construct mae columns names, and find corresponding param table cols */
    std::vector<std::string> maecols(nprops);
    std::vector<Id> propcols(nprops);
    for (int i=0; i<nprops; i++) {
        std::stringstream ss;
        ss << "ffio_c" << i+1;
        vdw.add_schema('r', ss.str());
        maecols[i]=ss.str();
        propcols[i] = table->propIndex(vdwprops[i]);
    }

    /* construct string keys for params */
    std::vector<std::string> vdwtypes;

    /* write the vdwtypes using its 1-based index as its name */
    for (unsigned i=0; i<table->paramCount(); i++) {
        std::stringstream ss;
        ss << i+1;
        std::string key = ss.str();
        vdwtypes.push_back(key);

        Destro& row = vdw.append();
        row["ffio_name"]=key;
        row["ffio_funct"]=funct;
        for (int j=0; j<nprops; j++) {
            row[maecols[j]] = table->paramTable()->value(i,propcols[j]).asFloat();
        }
    }

    /* now go back and updates ffio_sites */
    IdList terms = table->terms();
    for (unsigned i=0; i<terms.size(); i++) {
        Id id = terms[i];
        Id param = table->param(id);
        Id atom = table->atoms(id)[0];
        sites[atom+1]["ffio_vdwtype"] = vdwtypes[param];
    }
}

static void build_ff( SystemPtr mol, Destro& ffio_ff ) {
    std::vector<std::string> tables = mol->tableNames();
    for (unsigned i=0; i<tables.size(); i++) {
        const std::string& name = tables[i];
        //printf("processing %s\n", tables[i].c_str());
        TermTablePtr table = mol->table(name);
        if (name=="nonbonded") {
            build_nonbonded( mol, table, ffio_ff );
        } else { 
            build_tuple_table( mol, table, name, ffio_ff );
        }
    }
    /* special case for cmap tables */
    build_extra( mol, ffio_ff );
}

namespace desres { namespace msys {
    void ExportMAE( SystemPtr mol, String const& path, bool with_forcefield ) {

        Maeff M;
        std::ofstream out(path.c_str());

        /* create a single ct for the entire dms file */
        Destro& ct = M.new_block("f_m_ct");

        /* Get rid of any gaps in the atom list */
        if (mol->atomCount() != mol->maxAtomId()) mol=Clone(mol, mol->atoms());

        /* fill in the top-level ct stuff */
        build_ct_fields( mol, ct );

        /* add the atoms to the ct */
        //printf("processing atoms\n");
        build_m_atom( mol, ct );


        /* add the bonds to the ct */
        //printf("processing bonds\n");
        build_m_bond( mol, ct );

        if (with_forcefield){
            /* build sites and pseudos */
            Destro& ffio_ff = ct.new_block("ffio_ff");
            //printf("processing sites\n");
            build_sites( mol, ffio_ff );
            //printf("processing pseudos\n");
            build_pseudos( mol, ffio_ff );

            build_ff( mol, ffio_ff );
        }

        out << M;
    }
}}

