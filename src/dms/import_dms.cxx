#include "dms.hxx"
#include "../dms.hxx"
#include "../term_table.hxx"

#include <sstream>
#include <stdio.h>
#include <string.h>
#include <stdexcept>

#include <boost/algorithm/string.hpp> /* for boost::trim */
#include <boost/scoped_ptr.hpp>

using namespace desres::msys;

typedef std::set<String> KnownSet;

namespace {

    struct ResKey {
        int     resnum;
        String  resname;
        String  chain;

        ResKey() {}
        ResKey(int num, const String& name, const String& chn)
        : resnum(num), resname(name), chain(chn) {}

        bool operator<(const ResKey& r) const {
            if (resnum!=r.resnum) return resnum<r.resnum;
            int rc = strcmp(resname.c_str(), r.resname.c_str());
            if (rc) return rc<0;
            return strcmp(chain.c_str(), r.chain.c_str())<0;
        }
    };
}

static bool is_pN(const char * s) {
    if (s[0]!='p') return false;
    while (*(++s)) {
        if (*s<'0' || *s>'9') return false;
    }
    return true;
}

static void read(dms_reader_t* r, int col, ValueRef val) {
    switch (val.type()) {
        case FloatType: 
            val=dms_reader_get_double(r,col); 
            break;
        case IntType: 
            val=dms_reader_get_int(r,col); 
            break;
        default:
        case StringType: 
            val=dms_reader_get_string(r,col); 
            break;
    }
}

static IdList read_params( dms_reader_t* r, ParamTablePtr p,
                           bool ignore_ids = true ) {
    int i,n = dms_reader_column_count(r);
    int idcol=-1;
    for (i=0; i<n; i++) {
        const char * prop = dms_reader_column_name(r,i);
        /* ignore id, assuming param ids are 0-based */
        if (ignore_ids && !strcmp(prop, "id")) {
            idcol=i;
            continue;
        }
        p->addProp(prop, dms_reader_column_type(r,i));
    }
    IdList idmap;
    for (; r; dms_reader_next(&r)) {
        Id param = p->addParam();
        int j=0;
        for (i=0; i<n; i++) {
            if (i==idcol) {
                Id id = dms_reader_get_int(r,i);
                while (idmap.size()<id) idmap.push_back(BadId);
                idmap.push_back(param);
                continue;
            }
            read(r,i,p->value(param,j++));
        }
    }
    return idmap;
}

/* some people hack up their own dms file and put the terms and params into
 * one table instead of splitting into two plus a view.  Handle that case,
 * but refuse to deal with alchemical systems. */
static void read_view_table( dms_t* dms, 
                             const IdList& gidmap,
                             System& sys, 
                             const std::string& category,
                             const std::string& table,
                             KnownSet& known ) {

    dms_reader_t* r;

    /* find the terms table */
    known.insert(table);
    dms_fetch(dms,table.c_str(),&r);
    if (!r) return;

    std::vector<int> acols, pcols;
    ParamTablePtr params = ParamTable::create();
    int i,n = dms_reader_column_count(r);
    for (i=0; i<n; i++) {
        const char* prop = dms_reader_column_name(r,i);
        if (is_pN(prop)) acols.push_back(i);
        else {
            params->addProp(prop, dms_reader_column_type(r,i));
            pcols.push_back(i);
        }
    }

    TermTablePtr terms = sys.addTable(table, acols.size(), params);
    terms->category = category;
    IdList atoms(acols.size());
    for (; r; dms_reader_next(&r)) {
        for (unsigned i=0; i<acols.size(); i++) {
            atoms[i] = gidmap.at(dms_reader_get_int(r,acols[i]));
        }
        Id param = params->addParam();
        for (unsigned i=0; i<pcols.size(); i++) {
            read(r,pcols[i],params->value(param,i));
        }
        terms->addTerm(atoms, param);
    }
}

static void read_table( dms_t* dms, 
                        const IdList& gidmap,
                        System& sys, 
                        const std::string& category,
                        const std::string& table,
                        KnownSet& known ) {

    std::string term_table = table + "_term";
    std::string param_table = table + "_param";
    bool alchemical = false;
    if (table.substr(0,11)=="alchemical_") {
        param_table = param_table.substr(11);
        alchemical = true;
    }
    if (!(dms_has_table(dms, term_table.c_str()) &&
          dms_has_table(dms, param_table.c_str()))) {
        if (alchemical) {
            std::stringstream ss;
            ss << category << " table '" << table << "' is not well formatted.";
            throw std::runtime_error(ss.str());
        }
        return read_view_table(dms, gidmap, sys, category, table, known);
    }

    dms_reader_t* r;

    /* find the terms table */
    known.insert(term_table);
    known.insert(param_table);
    dms_fetch(dms,term_table.c_str(),&r);
    if (!r) return;

    /* get the number of particle columns */
    std::vector<int> cols;
    /* everything not an atom or param is an extra term property */
    typedef std::map<int,std::pair<String,ValueType> > ExtraMap;
    ExtraMap extra;
    int paramcol=-1, paramBcol=-1;
    int i,n = dms_reader_column_count(r);
    for (i=0; i<n; i++) {
        const char* prop = dms_reader_column_name(r,i);
        if (is_pN(prop)) cols.push_back(i);
        else if (!strcmp(prop, "param")) paramcol=i;
        else if (!strcmp(prop, "paramA")) paramcol=i;
        else if (!strcmp(prop, "paramB")) paramBcol=i;
        else {
            extra[i]=std::make_pair(
                    std::string(prop),dms_reader_column_type(r,i));
        }
    }
    unsigned natoms = cols.size();

    /* create and fill a TermTable entry in the System */
    TermTablePtr terms;
    if (alchemical) {
        terms = sys.addTable(table.substr(11), natoms);
    } else {
        terms = sys.addTable(table, natoms);
    }
    terms->category = category;
    for (ExtraMap::const_iterator i=extra.begin(); i!=extra.end(); ++i) {
        terms->addTermProp(i->second.first, i->second.second);
    }

    /* fill the ParamTable for the TermTable, but only if it hasn't already
     * been filled. */
    IdList idmap;
    {
        dms_reader_t* rp;
        dms_fetch(dms,param_table.c_str(),&rp);
        if (terms->params()->paramCount()) {
            /* must have already been filled by the alchemical or 
             * non-alchemical version of this table.  Make a dummy 
             * ParamTablePtr so that we can get the idmap. */
            idmap = read_params(rp, ParamTable::create());
        } else {
            idmap = read_params(rp, terms->params());
        }
    }
    //printf("created %s with %d params, %d props\n",
        //param_table.c_str(), params->paramCount(), params->propCount());

    IdList atoms(natoms);
    for (; r; dms_reader_next(&r)) {
        for (unsigned i=0; i<natoms; i++) {
            atoms[i] = gidmap.at(dms_reader_get_int(r,cols[i]));
        }
        Id param = paramcol>=0 ? dms_reader_get_int(r,paramcol) : BadId;
        if (!bad(param)) param=idmap.at(param);
        Id term = terms->addTerm(atoms, param);
        if (paramBcol>=0) {
            Id paramB = dms_reader_get_int(r,paramBcol);
            terms->setParamB(term,idmap.at(paramB));
        }
        /* extra properties */
        Id j=0;
        for (ExtraMap::const_iterator e=extra.begin(); e!=extra.end(); ++e) {
            read(r,e->first,terms->termPropValue(term, j++));
        }
    }
    //printf("created %s with %d atoms, %d extra properties, %d terms\n",
        //term_table.c_str(), 
        //terms->atomCount(), 
        //terms->props().propCount(),
        //terms->termCount());
}

static void read_metatables(dms_t* dms, const IdList& gidmap, System& sys,
                             KnownSet& known) {
    static const char * categories[] = { 
        "bond", "constraint", "virtual", "polar" 
    };
    dms_reader_t* r;
    for (unsigned i=0; i<sizeof(categories)/sizeof(categories[0]); i++) {
        std::string category = categories[i];
        std::string metatable = category + "_term";
        known.insert(metatable);
        if (dms_fetch(dms,metatable.c_str(),&r)) {
            int col=dms_reader_column(r,"name");
            for (; r; dms_reader_next(&r)) {
                std::string table = dms_reader_get_string(r,col);
                read_table( dms, gidmap, sys, category, table, known );
            }
        }
    }
}

static void 
read_nonbonded( dms_t* dms, System& sys, const std::vector<Id>& nbtypes,
                const std::map<Id,Id>& nbtypesB, KnownSet& known ) {

    known.insert("nonbonded_param");
    dms_reader_t* r;
    if (!dms_fetch(dms,"nonbonded_param",&r)) return;

    TermTablePtr terms = sys.addTable("nonbonded", 1);
    terms->category="nonbonded";
    IdList idmap = read_params(r, terms->params());

    IdList atoms(1);
    unsigned i,n = nbtypes.size();
    for (i=0; i<n; i++) {
        atoms[0]=i;
        Id nb = nbtypes[i];
        if (!bad(nb)) {
            try {
                nb = idmap.at(nb);
            }
            catch (std::exception& e) {
                std::stringstream ss;
                ss << "ImportDMS: particle " << i << " has invalid nbtype " << nb;

                throw std::runtime_error(ss.str());
            }
        }
        Id term = terms->addTerm(atoms, nb);
        std::map<Id,Id>::const_iterator iter = nbtypesB.find(atoms[0]);
        if (iter!=nbtypesB.end()) {
            terms->setParamB(term,iter->second);
        }
    }
}

static void
read_exclusions(dms_t* dms, const IdList& gidmap, System& sys, KnownSet& known) {
    known.insert("exclusion");
    /* some dms writers export exclusions with a term and param table. */
    known.insert("exclusion_term");
    known.insert("exclusion_param");

    dms_reader_t* r;
    if (!dms_fetch(dms,"exclusion",&r)) return;

    TermTablePtr terms = sys.addTable("exclusion", 2);
    terms->category="exclusion";
    IdList atoms(2);

    for (; r; dms_reader_next(&r)) {
        atoms[0] = gidmap.at(dms_reader_get_int(r,0));
        atoms[1] = gidmap.at(dms_reader_get_int(r,1));
        terms->addTerm(atoms,-1);
    }
}

static void
read_nbinfo(dms_t* dms, System& sys, KnownSet& known) {
    dms_reader_t *r;
    known.insert("nonbonded_info");
    if (dms_fetch(dms,"nonbonded_info",&r)) {
        int funct_col = dms_reader_column(r,"vdw_funct");
        int rule_col = dms_reader_column(r, "vdw_rule");
        if (funct_col>=0) {
            sys.nonbonded_info.vdw_funct = 
                dms_reader_get_string(r,funct_col);
        }
        if (rule_col>=0) {
            sys.nonbonded_info.vdw_rule =
                dms_reader_get_string(r,rule_col);
        }
    }
    dms_reader_free(r);
}

static void
read_cell(dms_t* dms, System& sys, KnownSet& known) {
    dms_reader_t* r;
    known.insert("global_cell");
    if (dms_fetch(dms,"global_cell",&r)) {
        int col[3];
        col[0]=dms_reader_column(r,"x");
        col[1]=dms_reader_column(r,"y");
        col[2]=dms_reader_column(r,"z");
        GlobalCell& cell = sys.global_cell;
        for (int i=0; i<3; i++) cell.A[i]=dms_reader_get_double(r,col[i]);
        dms_reader_next(&r);
        for (int i=0; i<3; i++) cell.B[i]=dms_reader_get_double(r,col[i]);
        dms_reader_next(&r);
        for (int i=0; i<3; i++) cell.C[i]=dms_reader_get_double(r,col[i]);
        dms_reader_next(&r);
        if (r) {
            throw std::runtime_error("global_cell table has too many rows");
        }
    }
}

static void read_extra( dms_t* dms, System& sys, const KnownSet& known) {
    dms_reader_t* r;
    if (dms_fetch(dms,"sqlite_master",&r)) {
        int NAME = dms_reader_column(r,"name");
        int TYPE = dms_reader_column(r,"type");
        if (NAME<0 || TYPE<0) {
            throw std::runtime_error("malformed sqlite_master table");
        }
        for (; r; dms_reader_next(&r)) {
            if (!strcmp(dms_reader_get_string(r,TYPE), "table")) {
                std::string extra = dms_reader_get_string(r,NAME);
                if (known.count(extra)) continue;
                dms_reader_t* p;
                if (dms_fetch(dms, extra.c_str(), &p)) {
                    ParamTablePtr ptr = ParamTable::create();
                    sys.addExtra(extra, ptr);
                    read_params(p, ptr, false);
                    //printf("extra %s with %d params, %d props\n",
                        //extra.c_str(), ptr->paramCount(), ptr->propCount());
                }
            }
        }
    }
}

static void read_provenance( dms_t* dms, System& sys, KnownSet& known) {
    known.insert("provenance");
    dms_reader_t* r;
    if (dms_fetch(dms, "provenance", &r)) {
        int timestamp = dms_reader_column(r, "timestamp");
        int user = dms_reader_column(r, "user");
        int workdir = dms_reader_column(r, "workdir");
        int cmdline = dms_reader_column(r, "cmdline");

        for (; r; dms_reader_next(&r)) {
            Provenance p;
            if (timestamp>=0) p.timestamp = dms_reader_get_string(r, timestamp);
            if (user>=0)      p.user      = dms_reader_get_string(r, user);
            if (workdir>=0)   p.workdir   = dms_reader_get_string(r, workdir);
            if (cmdline>=0)   p.cmdline   = dms_reader_get_string(r, cmdline);
            sys.addProvenance(p);
        }
    }
}

static void
read_alchemical_particle( dms_t* dms, System& sys, 
                          IdList& nbtypes, std::map<Id,Id>& nbtypesB,
                          KnownSet& known ) {
    known.insert("alchemical_particle");
    dms_reader_t* r;
    if (dms_fetch(dms,"alchemical_particle",&r)) {
        int P0 = dms_reader_column(r,"p0");
        int MOIETY = dms_reader_column(r,"moiety");
        int TYPEA = dms_reader_column(r,"nbtypeA");
        int TYPEB = dms_reader_column(r,"nbtypeB");
        int CHARGEA = dms_reader_column(r,"chargeA");
        int CHARGEB = dms_reader_column(r,"chargeB");

        if (P0<0 || MOIETY<0 || TYPEA<0 || TYPEB<0
                 || CHARGEA<0 || CHARGEB<0) {
            throw std::runtime_error("malformed alchemical_particle table");
        }
        for (; r; dms_reader_next(&r)) {
            Id id = dms_reader_get_int(r,P0);
            if (!sys.hasAtom(id)) {
                std::stringstream ss;
                ss << "alchemical_particle table has bad p0 '" << id << "'";
                throw std::runtime_error(ss.str());
            }
            atom_t& atm = sys.atom(id);
            atm.alchemical = true;
            atm.moiety = dms_reader_get_int(r,MOIETY);
            atm.charge = dms_reader_get_double(r,CHARGEA);
            atm.chargeB = dms_reader_get_double(r,CHARGEB);
            nbtypes.at(id) = dms_reader_get_int(r,TYPEA);
            nbtypesB[id] = dms_reader_get_int(r,TYPEB);
        }
    }
}
                                       
static SystemPtr import_dms( dms_t* dms, bool structure_only ) {

    SystemPtr h = System::create();
    System& sys = *h;
    dms_reader_t * r;

    Id chnid = BadId;
    Id resid = BadId;

    typedef std::map<String,Id> ChnMap;
    typedef std::map<ResKey,Id> ResMap;
    ResMap resmap;
    ChnMap chnmap;
    IdList nbtypes;
    IdList gidmap; /* map dms gids to msys ids */
    IdList ignored_gids;
    std::map<Id,Id> nbtypesB;
    
    dms_fetch(dms, "particle", &r);
    if (!r) {
        throw std::runtime_error("Missing particle table");
    }

    int CHAIN = dms_reader_column(r,"chain");
    int RESNAME = dms_reader_column(r,"resname");
    int RESID = dms_reader_column(r,"resid");
    int X = dms_reader_column(r,"x");
    int Y = dms_reader_column(r,"y");
    int Z = dms_reader_column(r,"z");
    int VX = dms_reader_column(r,"vx");
    int VY = dms_reader_column(r,"vy");
    int VZ = dms_reader_column(r,"vz");
    int MASS = dms_reader_column(r,"mass");
    int ANUM = dms_reader_column(r,"anum");
    int NAME = dms_reader_column(r,"name");
    int NBTYPE = dms_reader_column(r,"nbtype");
    int GID = dms_reader_column(r,"id");
    int CHARGE = dms_reader_column(r,"charge");
    
    /* the rest of the columns are extra atom properties */
    std::set<int> handled;
    handled.insert(CHAIN);
    handled.insert(RESNAME);
    handled.insert(RESID);
    handled.insert(X);
    handled.insert(Y);
    handled.insert(Z);
    handled.insert(VX);
    handled.insert(VY);
    handled.insert(VZ);
    handled.insert(MASS);
    handled.insert(ANUM);
    handled.insert(NAME);
    handled.insert(NBTYPE);
    handled.insert(GID);
    handled.insert(CHARGE);

    typedef std::map<int,ValueType> ExtraMap;
    ExtraMap extra;
    for (int i=0, n=dms_reader_column_count(r); i<n; i++) {
        if (handled.count(i)) continue;
        extra[i]=dms_reader_column_type(r,i);
        sys.addAtomProp(dms_reader_column_name(r,i), extra[i]);
    }

    /* read the particle table */
    for (; r; dms_reader_next(&r)) {
        Id gid = dms_reader_get_int(r, GID);
        int anum = dms_reader_get_int(r, ANUM);
        if (structure_only && anum<=0) {
            ignored_gids.push_back(gid);
            continue;
        }

        /* start a new chain if necessary */
        const char * chainname = dms_reader_get_string(r,CHAIN);
        ChnMap::const_iterator cit = chnmap.find(chainname);
        if (cit==chnmap.end()) {
            chnid = sys.addChain();
            sys.chain(chnid).name = chainname;
            resid = BadId;
            chnmap[chainname]=chnid;
        } else {
            chnid = cit->second;
        }

        /* start a new residue if necessary */
        const char * resname = dms_reader_get_string(r, RESNAME);
        int resnum = dms_reader_get_int(r, RESID);
        std::pair<ResMap::iterator,bool> p;
        p = resmap.insert(std::make_pair(ResKey(resnum,resname,chainname),
                    resid));
        if (p.second) {
            /* new resname/resnum in this chain, so start a new residue. */
            resid = sys.addResidue(chnid);
            residue_t& res = sys.residue(resid);
            p.first->second = resid;
            res.name = resname;
            res.num = resnum;

            boost::trim(res.name);
        } else {
            /* use existing residue */
            resid = p.first->second;
        }

        /* add the atom */
        Id atmid = sys.addAtom(resid);
        atom_t& atm = sys.atom(atmid);
        atm.x = dms_reader_get_double(r, X);
        atm.y = dms_reader_get_double(r, Y);
        atm.z = dms_reader_get_double(r, Z);
        atm.vx = dms_reader_get_double(r, VX);
        atm.vy = dms_reader_get_double(r, VY);
        atm.vz = dms_reader_get_double(r, VZ);
        atm.mass = dms_reader_get_double(r, MASS);
        atm.atomic_number = anum;
        atm.name = dms_reader_get_string(r, NAME);
        atm.gid = gid;
        atm.charge = dms_reader_get_double(r, CHARGE);
        while (gidmap.size()<atm.gid) gidmap.push_back(BadId);
        gidmap.push_back(atmid);

        boost::trim(atm.name);

        /* extra atom properties */
        Id propcol=0;
        for (ExtraMap::const_iterator iter=extra.begin();
                iter!=extra.end(); ++iter, ++propcol) {
            int col = iter->first;
            ValueType type = iter->second;
            ValueRef ref = sys.atomPropValue(atmid, propcol);
            if (type==IntType) 
                ref=dms_reader_get_int(r,col);
            else if (type==FloatType)
                ref=dms_reader_get_double(r,col);
            else
                ref=dms_reader_get_string(r,col);
        }
        nbtypes.push_back(NBTYPE>=0 ? dms_reader_get_int(r,NBTYPE) : BadId);
    }

    dms_fetch(dms, "bond", &r);
    if (r) {
        int p0=dms_reader_column(r,"p0");
        int p1=dms_reader_column(r,"p1");
        int o = dms_reader_column(r,"order");

        /* read the bond table */
        for (; r; dms_reader_next(&r)) {
            int ai = dms_reader_get_int(r, p0);
            int aj = dms_reader_get_int(r, p1);
            double order = o<0 ? 1 : dms_reader_get_double(r, o);
            Id id = BadId;
            if (structure_only && (
                        std::binary_search(
                            ignored_gids.begin(), ignored_gids.end(), ai) ||
                        std::binary_search(
                            ignored_gids.begin(), ignored_gids.end(), aj))) {
                continue;
            }
            try {
                Id from = gidmap.at(ai);
                Id to   = gidmap.at(aj);
                //printf("ai %d aj %d from %d to %d max id %d\n", 
                        //ai, aj, (int)from, (int)to, (int)sys.maxAtomId());
                id = sys.addBond(from,to);
            }
            catch (std::exception& e) {
                std::stringstream ss;
                ss << "Failed adding bond (" << ai << ", " << aj << ")"
                   << ": " << std::endl << e.what();
                throw std::runtime_error(ss.str());
            }
            sys.bond(id).order = order;
        }
    }

    sys.updateFragids();

    KnownSet known;
    known.insert("particle");
    known.insert("bond");

    read_cell(dms, sys, known);
    read_alchemical_particle( dms, sys, nbtypes, nbtypesB, known );
    read_provenance(dms, sys, known);

    if (!structure_only) {
        read_metatables(dms, gidmap, sys, known);
        read_nonbonded(dms, sys, nbtypes, nbtypesB, known);
        read_exclusions(dms, gidmap, sys, known);
        read_nbinfo(dms, sys, known);
        read_extra(dms, sys, known);
    }

    return h;
}

SystemPtr desres::msys::ImportDMS(const std::string& path, 
                                  bool structure_only) {
    dms_t * dms = NULL;
    SystemPtr sys;
    try {
        dms = dms_read(path.c_str());
        sys = import_dms(dms, structure_only);
    }
    catch (std::exception& e) {
        if (dms) dms_close(dms);
        std::stringstream ss;
        ss << "Error opening dms file at '" << path << "': " << e.what();
        throw std::runtime_error(ss.str());
    }
    dms_close(dms);
    return sys;
}

SystemPtr desres::msys::sqlite::ImportDMS(sqlite3* db,
                                  bool structure_only) {
    boost::scoped_ptr<dms_t> p(new dms_t(db));
    SystemPtr sys = import_dms(p.get(), structure_only);
    return sys;
}

