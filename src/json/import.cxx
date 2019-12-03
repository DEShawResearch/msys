#include "../json.hxx"
#include "../analyze.hxx"
#include <sys/stat.h>
#include <fcntl.h>

#if defined __has_include
#  if __has_include (<rapidjson/document.h>)
#    include <rapidjson/document.h>
#    include <rapidjson/reader.h>
#    define MSYS_WITH_RAPID_JSON
using namespace rapidjson;
#  endif
#endif


using namespace desres;

using msys::System;
using msys::SystemPtr;
using msys::ParamTablePtr;
using msys::ParamTable;
using msys::Id;

#if defined(MSYS_WITH_RAPID_JSON)

static std::shared_ptr<char> slurp(const char* path) {

    int fd=open(path, O_RDONLY);
    if (fd<0) {
        MSYS_FAIL("Reading file at '" << path << "': " << strerror(errno));
    }
    struct stat statbuf[1];
    if (fstat(fd, statbuf)!=0) {
        int _errno = errno;
        ::close(fd);
        MSYS_FAIL("Getting size of file at '" << path << "': " << strerror(_errno));
    }

    ssize_t tmpsize = statbuf->st_size;
    if (tmpsize==0) {
        close(fd);
        MSYS_FAIL("file at '" << path << "' has zero size");
    }
    char* tmpbuf = (char *)malloc(tmpsize+1);
    if (!tmpbuf) {
        close(fd);
        MSYS_FAIL("Failed to allocate read buffer for file at '" << path
           << "' of size " << tmpsize);
    }
    ssize_t sz = tmpsize;
    char* ptr = tmpbuf;
    while (sz) {
        errno = 0;
        ssize_t rc = ::read(fd, ptr, sz);
        if (rc<0 || (rc==0 && errno!=0)) {
            std::string errmsg = strerror(errno);
            close(fd);
            free(tmpbuf);
            MSYS_FAIL("Error reading file contents at " << path 
                    << ": " << errmsg);
        }
        sz -= rc;
        ptr += rc;
    }
    *ptr++ = '\0';
    close(fd);
    return std::shared_ptr<char>(tmpbuf, free);
}

static void read_chains(Document const& d, SystemPtr mol) {
    auto& names = d["names"];
    auto& chains = d["chains"];
    if (chains.ObjectEmpty()) {
        mol->addChain();
        return;
    }
    
    // required fields
    auto& name = chains["name"];

    // optional fields
    auto const& segid = chains.FindMember("segid");

    for (Id i=0, n=chains["count"].GetInt(); i<n; i++) {
        auto& chn = mol->chainFAST(mol->addChain());
        chn.name = names[name[i].GetInt()].GetString();
        if (segid != chains.MemberEnd()) {
            chn.segid = names[segid->value[i].GetInt()].GetString();
        }
    }
}

static void read_residues(Document const& d, SystemPtr mol) {
    auto& names = d["names"];
    auto& residues = d["residues"];
    if (residues.ObjectEmpty()) {
        mol->addResidue(0);
        return;
    }
    
    // required fields
    auto& chain = residues["chain"];
    auto& resid = residues["resid"];
    auto& name = residues["name"];

    // optional fields
    auto const& insertion = residues.FindMember("insertion");

    for (Id i=0, n=residues["count"].GetInt(); i<n; i++) {
        auto& res = mol->residueFAST(mol->addResidue(chain[i].GetInt()));
        res.name = names[name[i].GetInt()].GetString();
        res.resid = resid[i].GetInt();

        if (insertion != residues.MemberEnd()) {
            res.insertion = names[insertion->value[i].GetInt()].GetString();
        }
    }
}

static msys::ValueType parse_type(const char *s) {
    using namespace msys;
    switch (*s) {
        case 'i': return IntType;
        case 'f': return FloatType;
        default:;
        case 's': return StringType;
    }
}

static void read_tags(Value const& tags, ParamTablePtr params, Value const& names) {
    using msys::IntType;
    using msys::FloatType;
    using msys::StringType;
    for (auto& m : tags.GetObject()) {
        auto type = parse_type(m.value["type"].GetString());
        Id propid = params->addProp(m.name.GetString(), type);
        Value const& ids = m.value["ids"];
        Value const& vals = m.value["vals"];
        for (Id i=0, n=m.value["count"].GetInt(); i<n; i++) {
            Id id = ids[i].GetInt();
            while (params->paramCount() < id) params->addParam();
            auto ref = params->value(id, propid);
            switch (type) {
                case IntType: 
                    ref = vals[i].GetInt();
                    break;
                case FloatType:
                    ref = vals[i].GetDouble();
                    break;
                case StringType:
                    ref = names[vals[i].GetInt()].GetString();
                    break;
            }
        }
    }
}

static void read_cell(Document const& d, SystemPtr mol) {
    auto& cell = d["cell"];
    if (cell.Empty()) return;
    double* dst = mol->global_cell[0];
    for (int i=0; i<9; i++) {
        dst[i] = cell[i].GetDouble();
    }
}

static void read_particles(Document const& d, SystemPtr mol) {
    auto& names = d["names"];
    auto& particles = d["particles"];
    const Id natoms = particles["count"].GetInt();
    // required fields
    auto& anum = particles["atomic_number"];
    auto& fc = particles["formal_charge"];
    auto& name = particles["name"];

    // optional fields
    auto const& pos = particles.FindMember("position");
    auto const& vel = particles.FindMember("velocity");
    auto const& residue = particles.FindMember("residue");
    auto const& mass = particles.FindMember("mass");
    auto const& charge = particles.FindMember("charge");
    

    for (Id i=0; i<natoms; i++) {
        Id res = residue != particles.MemberEnd() ? residue->value[i].GetInt() : 0;
        auto& atm = mol->atomFAST(mol->addAtom(res));
        atm.name = names[name[i].GetInt()].GetString();
        if (pos != particles.MemberEnd()) {
            atm.x = pos->value[3*i  ].GetDouble();
            atm.y = pos->value[3*i+1].GetDouble();
            atm.z = pos->value[3*i+2].GetDouble();
        }
        if (vel != particles.MemberEnd()) {
            atm.vx = vel->value[3*i  ].GetDouble();
            atm.vy = vel->value[3*i+1].GetDouble();
            atm.vz = vel->value[3*i+2].GetDouble();
        }
        atm.atomic_number = anum[i].GetInt();
        atm.formal_charge = fc[i].GetInt();

        if (mass != particles.MemberEnd()) atm.mass = mass->value[i].GetDouble();
        if (charge != particles.MemberEnd()) atm.charge = charge->value[i].GetDouble();
    }
    auto const& tags = particles.FindMember("tags");
    if (tags != particles.MemberEnd()) {
        read_tags(tags->value, mol->atomProps(), names);
    }

}

static void read_bonds(Document const& d, SystemPtr mol) {
    auto& bonds = d["bonds"];
    auto& p = bonds["particles"];
    auto& order = bonds["order"];
    for (Id i=0, n=bonds["count"].GetInt(); i<n; i++) {
        auto& bond = mol->bondFAST(mol->addBond(
            p[2*i].GetInt(),
            p[2*i+1].GetInt()));
        bond.order = order[i].GetInt();
    }
}

static void read_params(Value const& val, Value const& names, ParamTablePtr params) {
    const Id nparams = val["count"].GetInt();
    for (Id i=0; i<nparams; i++) params->addParam();
    for (auto const& m : val["props"].GetObject()) {
        auto type = parse_type(m.value["type"].GetString());
        auto vals = m.value["vals"].GetArray();
        Id j=params->addProp(m.name.GetString(), type);
        switch (type) {
            case msys::IntType:
                for (Id i=0; i<nparams; i++) {
                    params->value(i,j) = vals[i].GetInt();
                }
                break;
            case msys::FloatType:
                for (Id i=0; i<nparams; i++) {
                    params->value(i,j) = vals[i].GetDouble();
                }
                break;
            case msys::StringType:
                for (Id i=0; i<nparams; i++) {
                    params->value(i,j) = names[vals[i].GetInt()].GetString();
                }
                break;
        }
    }
}

static void read_aux(Document const& d, SystemPtr mol) {
    auto& names = d["names"];
    auto const& tables = d.FindMember("aux");
    if (tables == d.MemberEnd()) return;
    for (auto& m : tables->value.GetObject()) {
        auto params = ParamTable::create();
        read_params(m.value, names, params);
        mol->addAuxTable(m.name.GetString(), params);
    }
}

static void read_tables(Document const& d, SystemPtr mol) {
    auto& names = d["names"];
    auto const& tables = d.FindMember("tables");
    if (tables == d.MemberEnd()) return;
    for (auto& m : tables->value.GetObject()) {
        auto table = mol->addTable(m.name.GetString(), m.value["arity"].GetInt());
        auto attrs = m.value.FindMember("attrs");
        if (attrs != m.value.MemberEnd()) {
            auto const& category = attrs->value.FindMember("category");
            if (category != attrs->value.MemberEnd()) {
                table->category = msys::parse(category->value.GetString());
            }
            if (table->name() == "nonbonded") {
                auto const& rule = attrs->value.FindMember("vdw_rule");
                if (rule != attrs->value.MemberEnd()) {
                    mol->nonbonded_info.vdw_rule = rule->value.GetString();
                    mol->nonbonded_info.vdw_funct = "vdw_12_6";
                }
            }
        }
        read_params(m.value["param"], names, table->params());
        auto& terms = m.value["terms"];
        auto& particles = terms["particles"];
        auto& params = terms["params"];
        msys::IdList atoms(table->atomCount());
        Id particle_index = 0;
        for (Id i=0, n=terms["count"].GetInt(); i<n; i++) {
            for (Id j=0, m=table->atomCount(); j<m; j++) {
                atoms[j] = particles[particle_index++].GetInt();
            }
            table->addTerm(atoms, params[i].GetInt());
        }
        auto& tags = m.value["tags"];
        read_tags(tags, table->props(), names);
    }
}


static SystemPtr import_json(Document const& d) {
    auto mol = System::create();
    read_chains(d, mol);
    read_residues(d, mol);
    read_particles(d, mol);
    read_cell(d, mol);
    read_bonds(d, mol);
    read_tables(d, mol);
    read_aux(d, mol);
    msys::Analyze(mol);
    return mol;
}

namespace desres { namespace msys {

    SystemPtr ImportJson(std::string const& path) {
        Document d;
        auto json = slurp(path.data());
        d.Parse<kParseFullPrecisionFlag>(json.get());
        SystemPtr mol = import_json(d);
        return mol;
    }

    SystemPtr ParseJson(const char* text) {
        Document d;
        d.Parse<kParseFullPrecisionFlag>(text);
        SystemPtr mol = import_json(d);
        return mol;
    }

}}

#else
#warning "rapidjson not available; no json import support"
namespace desres { namespace msys {

    SystemPtr ImportJson(std::string const& path) {
        MSYS_FAIL("msys compiled without rapidjson support; json import not available");
        return nullptr;
    }

    SystemPtr ParseJson(const char* text) {
        MSYS_FAIL("msys compiled without rapidjson support; json import not available");
        return nullptr;
    }

}}

#endif