#include "dms.hxx"
#include "../dms.hxx"
#include "../term_table.hxx"

#include <sstream>
#include <stdio.h>
#include <string.h>
#include <stdexcept>
#include <boost/scoped_ptr.hpp>

using namespace desres::msys;

static const char* str(const ValueType& t) {
    return t==IntType   ? "integer" :
           t==FloatType ? "float" :
                          "string";
}

static void write(const ValueRef& ref, int col, dms_writer_t* w) {
    switch (ref.type()) {
        case IntType: 
            dms_writer_bind_int(w,col,ref.asInt()); 
            break;
        case FloatType: 
            dms_writer_bind_double(w,col,ref.asFloat()); 
            break;
        default:
        case StringType:
            dms_writer_bind_string(w,col,ref.asString().c_str());
            break;
    }
}

typedef std::pair<Id,Id> NbType;
typedef std::map<Id,NbType> NbMap;
NbMap fetch_nbtypes(const System& sys) {
    TermTablePtr table;
    std::vector<String> tables = sys.tableNames();
    for (unsigned i=0; i<tables.size(); i++) {
        TermTablePtr t = sys.table(tables[i]);
        if (t->category=="nonbonded") {
            if (table) {
                throw std::runtime_error("multiple nonbonded tables found");
            }
            if (t->atomCount()!=1) {
                throw std::runtime_error("nonbonded table has atomCount!=1");
            }
            table=t;
        }
    }
    NbMap nbtypes;
    if (table) {
        IdList ids=table->terms();
        for (Id i=0, n=ids.size(); i<n; i++) {
            Id termid=ids[i];
            Id atomid=table->atoms(termid)[0];
            Id type  = table->param(termid);
            Id typeB = table->paramB(termid);
            if (bad(typeB)) typeB=type;
            nbtypes[atomid].first =type;
            nbtypes[atomid].second=typeB;
        }
    }
    return nbtypes;
}

static void export_alchemical_particles(const System& sys, 
                                        const IdList& alchemical_ids, 
                                        NbMap& nbtypes,
                                        dms_t* dms) {
    if (!alchemical_ids.size()) return;
    dms_exec(dms, "create table alchemical_particle (\n"
                  "  p0 integer primary key,\n"
                  "  moiety integer, \n"
                  "  nbtypeA integer, \n" 
                  "  nbtypeB integer, \n"
                  "  chargeA float,\n"
                  "  chargeB float)");
    dms_writer_t* w;
    dms_insert(dms, "alchemical_particle", &w);
    for (Id i=0; i<alchemical_ids.size(); i++) {
        Id id = alchemical_ids[i];
        const atom_t& atom = sys.atom(id);
        dms_writer_bind_int(w,0,id);
        dms_writer_bind_int(w,1,atom.moiety);
        dms_writer_bind_int(w,2,nbtypes[id].first);
        dms_writer_bind_int(w,3,nbtypes[id].second);
        dms_writer_bind_double(w,4,atom.charge);
        dms_writer_bind_double(w,5,atom.chargeB);
        dms_writer_next(w);
    }
    dms_writer_free(w);
}

static void export_particles(const System& sys, const IdList& map, dms_t* dms) {

    NbMap nbtypes = fetch_nbtypes(sys);
    IdList ids = sys.atoms();
    IdList alchemical_ids;
    for (Id i=0; i<ids.size(); i++) {
        if (sys.atom(ids[i]).alchemical) {
            alchemical_ids.push_back(ids[i]);
        }
    }

    std::string sql = 
        "create table particle (\n"
        "  id integer primary key,\n"
        "  anum integer,\n"
        "  name text,\n" 
        "  x float,\n"
        "  y float,\n"
        "  z float,\n"
        "  vx float,\n"
        "  vy float,\n"
        "  vz float,\n"
        "  resname text,\n"
        "  resid integer,\n"
        "  chain text,\n"
        "  mass float,\n"
        "  charge float,\n";

    const Id nprops = sys.atomPropCount();
    for (Id i=0; i<nprops; i++) {
        sql += "  '" + sys.atomPropName(i) + "' ";
        sql += str(sys.atomPropType(i));
        sql += ",\n";
    }
    sql += "  nbtype integer\n);";
    dms_exec(dms, sql.c_str());

    dms_writer_t* w;
    dms_insert(dms,"particle", &w);
    dms_exec(dms, "begin");
    for (Id i=0, n=ids.size(); i<n; i++) {
        Id atm = ids[i];
        const atom_t& atom = sys.atom(atm);
        Id res = atom.residue;
        const residue_t& residue = sys.residue(res);
        Id chn = residue.chain;
        const chain_t& chain = sys.chain(chn);

        dms_writer_bind_int(w,    0, map[i]);
        dms_writer_bind_int(w,    1, atom.atomic_number);
        dms_writer_bind_string(w, 2, atom.name.c_str());
        dms_writer_bind_double(w, 3, atom.x);
        dms_writer_bind_double(w, 4, atom.y);
        dms_writer_bind_double(w, 5, atom.z);
        dms_writer_bind_double(w, 6, atom.vx);
        dms_writer_bind_double(w, 7, atom.vy);
        dms_writer_bind_double(w, 8, atom.vz);
        dms_writer_bind_string(w, 9, residue.name.c_str());
        dms_writer_bind_int   (w,10, residue.num);
        dms_writer_bind_string(w,11, chain.name.c_str());
        dms_writer_bind_double(w,12, atom.mass);
        dms_writer_bind_double(w,13, atom.charge);
        for (Id j=0; j<nprops; j++) {
            int col=14+j;

            /* *sigh* - the ParamTable::value() method is non-const,
             * and I don't feel like making a const version; thus this
             * hack. */
            ValueRef ref = const_cast<System&>(sys).atomPropValue(atm,j);
            write(ref, col, w);
        }
        NbMap::const_iterator nbiter=nbtypes.find(atm);
        if (nbiter!=nbtypes.end()) {
            dms_writer_bind_int(w,14+nprops,nbiter->second.first);
        } else {
            dms_writer_bind_null(w,14+nprops);
        }
        dms_writer_next(w);
    }
    dms_writer_free(w);
    export_alchemical_particles(sys, alchemical_ids, nbtypes, dms);
    dms_exec(dms, "commit");
}

static void export_bonds(const System& sys, const IdList& map, dms_t* dms) {
    std::string sql =
        "create table bond (\n"
        "  p0 integer,\n"
        "  p1 integer,\n"
        "  'order' float\n"
        ");";
    dms_exec(dms, sql.c_str());

    dms_writer_t* w;
    dms_insert(dms,"bond", &w);
    IdList ids = sys.bonds();
    dms_exec(dms, "begin");
    for (Id i=0, n=ids.size(); i<n; i++) {
        Id bnd = ids[i];
        const bond_t& bond = sys.bond(bnd);
        dms_writer_bind_int(w,0,map[bond.i]);
        dms_writer_bind_int(w,1,map[bond.j]);
        dms_writer_bind_double(w,2,bond.order);
        dms_writer_next(w);
    }
    dms_writer_free(w);
    dms_exec(dms, "commit");
}

static void export_terms(TermTablePtr table, const IdList& map, 
                         const std::string& tablename, dms_t* dms) {

    const Id natoms = table->atomCount();
    const Id nprops = table->termPropCount();
    std::string alcname("alchemical_"); alcname += tablename;
    std::stringstream ss;
    ss << "create table " << tablename << " (";
    for (Id i=0; i<natoms; i++) ss << "p" << i << " integer, ";
    for (Id i=0; i<nprops; i++) {
        ss << "'" << table->termPropName(i) << "' " 
           << str(table->termPropType(i)) << ", ";
    }
    ss << "param integer)";
    dms_exec(dms, ss.str().c_str());

    if (table->alchemical()) {
        std::stringstream ss;
        ss << "create table " << alcname << " (";
        for (Id i=0; i<natoms; i++) ss << "p" << i << " integer, ";
        for (Id i=0; i<nprops; i++) {
            ss << "'" << table->termPropName(i) << "' " 
               << str(table->termPropType(i)) << ", ";
        }
        ss << "paramA integer, paramB integer)";
        dms_exec(dms, ss.str().c_str());
    }

    dms_writer_t* regw, *alcw=NULL;
    dms_insert(dms,tablename.c_str(),&regw);
    if (table->alchemical()) {
        dms_insert(dms,alcname.c_str(),&alcw);
    }

    dms_exec(dms,"begin");
    IdList ids = table->terms();
    for (Id i=0,n=ids.size(); i<n; i++) {
        Id id=ids[i];
        Id paramB = table->paramB(id);

        /* see if this term is alchemical or not */
        dms_writer_t* w = regw;
        if (!bad(paramB)) w = alcw;

        /* write atom columns */
        IdList atoms = table->atoms(id);
        for (Id j=0; j<natoms; j++) dms_writer_bind_int(w,j,map[atoms[j]]);
        /* write extra atom properties */
        for (Id j=0; j<nprops; j++) {
            ValueRef val = table->termPropValue(id, j);
            write(val, j+natoms, w);
        }
        /* write param column, if there is a param */
        Id param = table->param(id);
        if (bad(param)) {
            dms_writer_bind_null(w,nprops+natoms);
        } else {
            dms_writer_bind_int(w,nprops+natoms,param);
        }
        if (w==alcw) dms_writer_bind_int(w,nprops+natoms+1,paramB);
        dms_writer_next(w);
    }
    if (regw) dms_writer_free(regw);
    if (alcw) dms_writer_free(alcw);
    dms_exec(dms,"commit");
}

static void export_params(ParamTablePtr params, const std::string& tablename,
                          dms_t* dms, bool with_id=true) {

    std::stringstream ss;
    const Id nprops=params->propCount();
    if (nprops==0 && !with_id) {
      /* no columns have been specified, so there are no relations and thus
       * nothing to do.  */
      return;
    }
    ss << "create table " << tablename << " (";
    for (Id i=0; i<nprops; i++) {
        ss << params->propName(i) << " " << str(params->propType(i));
        if (i!=nprops-1) ss << ", ";
    }
    if (with_id) {
        if (nprops>0) ss << ", ";
        ss << "id integer primary key)";
    } else {
        ss << ")";
    }

    dms_exec(dms, ss.str().c_str());
    dms_writer_t* w;
    dms_insert(dms,tablename.c_str(),&w);
    dms_exec(dms, "begin");
    for (Id i=0, n=params->paramCount(); i<n; i++) {
        for (Id j=0; j<nprops; j++) {
            ValueRef val = params->value(i,j);
            write(val, j, w);
        }
        if (with_id) {
            dms_writer_bind_int(w,nprops,i);
        }
        dms_writer_next(w);
    }
    dms_exec(dms, "commit");
    dms_writer_free(w);
}

static void export_view(TermTablePtr table, const std::string& name, dms_t *dms) {
    std::string termname = name + "_term";
    std::string paramname = name + "_param";
    std::vector<String> props, tprops;
    for (Id i=0; i<table->params()->propCount(); i++) {
        const std::string& prop = table->params()->propName(i);
        if (prop!="id") props.push_back(prop);
    }
    for (Id i=0; i<table->termPropCount(); i++) {
        const std::string& prop = table->termPropName(i);
        if (prop!="param") tprops.push_back(prop);
    }
    std::stringstream ss;
    ss << "create view " << name << " as \n" << "  select ";
    for (Id i=0; i<table->atomCount(); i++) {
        ss << "p" << (char)('0'+i);
        if (props.size() || tprops.size() || i!=table->atomCount()-1) {
            ss << ", ";
        }
    }
    ss << "\n";
    for (Id i=0; i<props.size(); i++) {
        ss << props[i];
        if (tprops.size() || i!=props.size()-1) ss << ", ";
    }
    for (Id i=0; i<tprops.size(); i++) {
        ss << tprops[i];
        if (i!=tprops.size()-1) ss << ", ";
    }
    ss << "  from " << paramname << "\n"
       << "  join " << termname << "\n"
       << "  on param=id";
    dms_exec(dms,ss.str().c_str());

    if (table->alchemical()) {
        std::stringstream ss;
        ss << "create view alchemical_" << name << " as \n" << "  select ";
        for (Id i=0; i<table->atomCount(); i++) {
            ss << "p" << (char)('0'+i) << ", ";
        }
        ss << "moiety, \n";
        for (Id i=0; i<props.size(); i++) {
            ss << "  A." << props[i] << " as " << props[i] << "A, ";
            ss << "B." << props[i] << " as " << props[i] << "B";
            if (i!=props.size()-1) ss << ", ";
            ss << "\n";
        }
        ss << "  from alchemical_" << termname << "\n"
           << "  join " << paramname << " as A on paramA=A.id\n"
           << "  join " << paramname << " as B on paramB=B.id\n";
        dms_exec(dms,ss.str().c_str());
    }
}

static void export_exclusion(TermTablePtr table, const IdList& map, dms_t* dms) {
    if (table->atomCount()!=2) {
        throw std::runtime_error("table with category exclusion has atomCount!=2");
    }
    dms_exec(dms, "create table exclusion (p0 integer, p1 integer)");
    IdList ids = table->terms();
    dms_writer_t* w;
    dms_insert(dms,"exclusion",&w);
    dms_exec(dms, "begin");
    for (Id i=0, n=ids.size(); i<n; i++) {
        IdList atoms = table->atoms(ids[i]);
        dms_writer_bind_int(w,0,map[atoms[0]]);
        dms_writer_bind_int(w,1,map[atoms[1]]);
        dms_writer_next(w);
    }
    dms_writer_free(w);
    dms_exec(dms, "commit");
}

static void export_nonbonded( TermTablePtr table, dms_t* dms) {
    if (table->atomCount()!=1) {
        throw std::runtime_error("table with category nonbonded has atomCount!=1");
    }
    export_params(table->params(), "nonbonded_param", dms);
}

static void export_meta( TermTablePtr table, const std::string& name, 
        dms_t* dms) {
    std::string sql("insert into ");
    sql += table->category;
    sql += "_term values ('";
    sql += name;
    sql += "')";
    dms_exec(dms, sql.c_str());
    if (table->alchemical()) {
        std::string sql("insert into ");
        sql += table->category;
        sql += "_term values ('alchemical_";
        sql += name;
        sql += "')";
        dms_exec(dms, sql.c_str());
    }
}

static void export_tables( const System& sys, const IdList& map, dms_t* dms) {
    dms_exec(dms, "create table bond_term (name text)");
    dms_exec(dms, "create table constraint_term (name text)");
    dms_exec(dms, "create table virtual_term (name text)");
    dms_exec(dms, "create table polar_term (name text)");
    std::vector<String> tables = sys.tableNames();
    for (unsigned i=0; i<tables.size(); i++) {
        const std::string& name = tables[i];
        TermTablePtr table = sys.table(name);
        if (table->category.size()==0) {
            std::stringstream ss;
            ss << "cannot export table '" << tables[i] << "' with no category";
            throw std::runtime_error(ss.str());

        } else if (table->category=="exclusion") {
            export_exclusion(table, map, dms);
        } else if (table->category=="nonbonded") {
            export_nonbonded(table, dms);
        } else {
            export_terms(table, map, name+"_term", dms);
            export_params(table->params(), name+"_param", dms);
            export_view(table, name, dms);
            export_meta(table, name, dms);
        }
    }
}

static void export_extra(const System& sys, dms_t* dms) {
    std::vector<String> extras = sys.extraNames();
    for (unsigned i=0; i<extras.size(); i++) {
        const std::string& name = extras[i];
        export_params(sys.extra(name), name, dms, false);
    }
}

static void export_nbinfo(const System& sys, dms_t* dms) {
    dms_exec(dms,"create table nonbonded_info (vdw_funct text, vdw_rule text)");
    dms_writer_t* w;
    dms_insert(dms,"nonbonded_info", &w);
    dms_writer_bind_string(w,0,sys.nonbonded_info.vdw_funct.c_str());
    dms_writer_bind_string(w,1,sys.nonbonded_info.vdw_rule.c_str());
    dms_writer_next(w);
    dms_writer_free(w);
}

static void export_cell(const System& sys, dms_t* dms) {
    const GlobalCell& cell = sys.global_cell;
    dms_exec(dms, 
            "create table global_cell (\n"
            "  id integer primary key,\n"
            "  x float, y float, z float)");
    dms_writer_t* w;
    dms_insert(dms,"global_cell",&w);
    dms_exec(dms, "begin");
    for (int i=0; i<3; i++) dms_writer_bind_double(w,i+1,cell.A[i]);
    dms_writer_next(w);
    for (int i=0; i<3; i++) dms_writer_bind_double(w,i+1,cell.B[i]);
    dms_writer_next(w);
    for (int i=0; i<3; i++) dms_writer_bind_double(w,i+1,cell.C[i]);
    dms_writer_next(w);
    dms_writer_free(w);
    dms_exec(dms, "commit");
}

static IdList map_gids(System const& sys) {
    IdList gids, ids = sys.atoms();
    IdList idmap(sys.maxAtomId());
    for (Id i=0, n=ids.size(); i<n; i++) {
        Id id = ids[i];
        Id gid = sys.atom(id).gid;
        idmap[id] = gid;
        gids.push_back(gid);
    }
    /* Ensure that the set of gids has no repeats. */
    std::sort(gids.begin(), gids.end());
    if (std::unique(gids.begin(), gids.end()) != gids.end()) {
        throw std::runtime_error("atom gids are not unique");
    }
    return idmap;
}

static void export_dms(SystemPtr h, dms_t* dms) {
    System& sys = *h;
    IdList atomidmap = map_gids(sys);
    export_particles(sys, atomidmap, dms);
    export_bonds(    sys, atomidmap, dms);
    export_tables(   sys, atomidmap, dms);
    export_extra(    sys,            dms);
    export_nbinfo(   sys,            dms);
    export_cell(     sys,            dms);
}

void desres::msys::ExportDMS(SystemPtr h, const std::string& path) {
    unlink(path.c_str());
    dms_t* dms = NULL;
    try {
        dms = dms_write(path.c_str());
        export_dms(h, dms);
    }
    catch(std::exception& e) {
        if (dms) dms_close(dms);
        std::stringstream ss;
        ss << "Error writing dms file at '" << path << "': " << e.what();
        throw std::runtime_error(ss.str());
    }
    dms_close(dms);
}

void desres::msys::sqlite::ExportDMS(SystemPtr h, sqlite3* db) {
    boost::scoped_ptr<dms_t> p(new dms_t(db));
    export_dms(h, p.get());
}

