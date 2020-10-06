#include "../json.hxx"
#include "../analyze.hxx"
#include <sys/stat.h>
#include <fcntl.h>

#if defined __has_include
#  if __has_include (<rapidjson/document.h>)
#    include <rapidjson/document.h>
#    include <rapidjson/reader.h>
#    include <rapidjson/error/en.h>
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

template<typename T>
void check_sizes(const char *fname, T const &a1, const char *name1, T const &a2, const char *name2) {
    if (a1.Size() != a2.Size()) {
        MSYS_FAIL("function " << fname << " expected " << name1 << " and " << name2 << " arrays to be the same size, found " << a1.Size() << " and " << a2.Size());
    }
}

template<typename T, typename S>
void check_size(const char *fname, T const &a1, const char *name1, S s) {
    if (a1.Size() != s) {
        MSYS_FAIL("function " << fname << " expected " << name1 << " arrays to be the size " << s << " found " << a1.Size());
    }
}

#define CHECK_SIZES(a1, name1, a2, name2) check_sizes(__FUNCTION__, a1, name1, a2, name2)
#define CHECK_SIZE(a1, name1, s) check_size(__FUNCTION__, a1, name1, s)

static const char* get_name(Value const& names, int64_t nameid) {
    if (!names.IsArray()) {
        if (nameid != 0) {
            MSYS_FAIL("Unable to translate non-zero nameid to string without 'names' array");
        }
        return "";
    }
    return names[nameid].GetString();
}

static void read_chains(Document const& d, SystemPtr mol) {
    auto o = d.GetObject();
    auto const& m = o.FindMember("chains");
    if (m == o.MemberEnd() || m->value.ObjectEmpty()) {
        mol->addChain();
        return;
    }
    auto& chains = m->value;

    // required fields
    auto& name = chains["name"];

    // optional fields
    auto const& segid = chains.FindMember("segid");
    if (segid != chains.MemberEnd()) CHECK_SIZES(name, "name", segid->value, "segid");

    auto& names = d["names"];
    for (Id i=0, n=name.Size(); i<n; i++) {
        auto& chn = mol->chainFAST(mol->addChain());
        chn.name = get_name(names, name[i].GetInt());
        if (segid != chains.MemberEnd()) {
            chn.segid = get_name(names, segid->value[i].GetInt());
        }
    }
}

static void read_residues(Document const& d, SystemPtr mol) {
    auto const &m = d.FindMember("residues");
    if (m == d.MemberEnd() || m->value.ObjectEmpty()) {
        mol->addResidue(0);
        return;
    }
    auto& residues = m->value;

    // required fields
    auto& chain = residues["chain"];
    auto& resid = residues["resid"];
    CHECK_SIZES(chain, "chain", resid, "resid");
    auto& name = residues["name"];
    CHECK_SIZES(chain, "chain", resid, "name");

    // optional fields
    auto const& insertion = residues.FindMember("insertion");
    if (insertion != residues.MemberEnd()) CHECK_SIZES(chain, "chain", insertion->value, "insertion");

    auto& names = d["names"];
    for (Id i=0, n=chain.Size(); i<n; i++) {
        auto& res = mol->residueFAST(mol->addResidue(chain[i].GetInt()));
        res.name = get_name(names, name[i].GetInt());
        res.resid = resid[i].GetInt();

        if (insertion != residues.MemberEnd()) {
            res.insertion = get_name(names, insertion->value[i].GetInt());
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
        auto type = parse_type(m.value["t"].GetString());
        Id propid = params->addProp(m.name.GetString(), type);
        Value const& ids = m.value["i"];
        Value const& vals = m.value["v"];
        CHECK_SIZES(ids, "i", vals, "v");
        for (Id i=0, n=ids.Size(); i<n; i++) {
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
                    ref = get_name(names, vals[i].GetInt());
                    break;
            }
        }
    }
}

static void read_cell(Document const& d, SystemPtr mol) {
    auto const &m = d.FindMember("cell");
    if (m == d.MemberEnd() || m->value.Empty()) {
        return;
    }
    auto& cell = m->value;

    double* dst = mol->global_cell[0];
    for (int i=0; i<9; i++) {
        dst[i] = cell[i].GetDouble();
    }
}

static void read_particles(Document const& d, SystemPtr mol) {
    auto const &names = d.FindMember("names");

    auto& particles = d["i"];
    Id natoms = msys::BadId;

    // optional fields
    auto const& anum = particles.FindMember("atomic_number");
    if (anum != particles.MemberEnd()) {
        natoms = anum->value.Size();
    }

    auto const& name = particles.FindMember("name");
    if (name != particles.MemberEnd()) {
        if (natoms == msys::BadId) natoms = name->value.Size();
        CHECK_SIZE(name->value, "name", natoms);
    }

    auto const& fc = particles.FindMember("formal_charge");
    if (fc != particles.MemberEnd()) {
        if (natoms == msys::BadId) natoms = fc->value.Size();
        CHECK_SIZE(fc->value, "formal_charge", natoms);
    }

    auto const& pos = particles.FindMember("position");
    if (pos != particles.MemberEnd()) {
        if (natoms == msys::BadId) natoms = pos->value.Size() / 3;
        CHECK_SIZE(pos->value, "position", natoms * 3);
    }

    auto const& vel = particles.FindMember("velocity");
    if (vel != particles.MemberEnd()) {
        if (natoms == msys::BadId) natoms = vel->value.Size() / 3;
        CHECK_SIZE(vel->value, "velocity", natoms * 3);
    }

    auto const& residue = particles.FindMember("residue");
    if (residue != particles.MemberEnd()) {
        if (natoms == msys::BadId) natoms = residue->value.Size();
        CHECK_SIZE(residue->value, "residue", natoms);
    }

    auto const& mass = particles.FindMember("m");
    if (mass != particles.MemberEnd()) {
        if (natoms == msys::BadId) natoms = mass->value.Size();
        CHECK_SIZE(mass->value, "m", natoms);
    }

    auto const& charge = particles.FindMember("c");
    if (charge != particles.MemberEnd()) {
        if (natoms == msys::BadId) natoms = charge->value.Size();
        CHECK_SIZE(charge->value, "c", natoms);
    }

    if (natoms == msys::BadId)
        MSYS_FAIL("Unable to find any atoms in the particles object!");

    for (Id i=0; i<natoms; i++) {
        Id res = residue != particles.MemberEnd() ? residue->value[i].GetInt() : 0;
        auto& atm = mol->atomFAST(mol->addAtom(res));
        if (name != particles.MemberEnd()) {
            atm.name = get_name(names->value, name->value[i].GetInt());
        }

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

        if (anum != particles.MemberEnd()) atm.atomic_number = anum->value[i].GetInt();
        if (fc != particles.MemberEnd()) atm.formal_charge = fc->value[i].GetInt();
        if (mass != particles.MemberEnd()) atm.mass = mass->value[i].GetDouble();
        if (charge != particles.MemberEnd()) atm.charge = charge->value[i].GetDouble();
    }
    auto const& tags = particles.FindMember("tags");
    if (tags != particles.MemberEnd()) {
        read_tags(tags->value, mol->atomProps(), names->value);
    }

}

static void read_bonds(Document const& d, SystemPtr mol) {
    auto& bonds = d["b"];
    auto& p = bonds["i"];
    Id nbonds = p.Size() / 2;

    auto const& order = bonds.FindMember("order");
    if (order != bonds.MemberEnd()) CHECK_SIZE(order->value, "order", nbonds);
    
    for (Id i=0, n=nbonds; i<n; i++) {
        auto& bond = mol->bondFAST(mol->addBond(
            p[2*i].GetInt(),
            p[2*i+1].GetInt()));

        if (order != bonds.MemberEnd()) {
            bond.order = order->value[i].GetInt();
        }
    }
}

static void read_params(Value const& val, Value const& names, ParamTablePtr params, const std::string &name) {
    Id nparams = msys::BadId;

    // search for the number of parameters
    auto const& count = val.FindMember("c");
    if (count != val.MemberEnd()) {
        nparams = count->value.GetInt();
    } else {
        for (auto const& m : val["p"].GetObject()) {
            auto const& valm = m.value.FindMember("v");
            if (valm == m.value.MemberEnd()) continue; // no vals, skip along

            nparams = m.value["v"].GetArray().Size(); // use the first one we find
            break;
        }
    }
    if (nparams == msys::BadId) {
        MSYS_FAIL("Failed to determine number of parameters for " << name);
    }
    // instantiate the parameters
    for (Id i=0; i<nparams; i++)  params->addParam();

    // now populate them with the appropriate values
    for (auto const& m : val["p"].GetObject()) {
        auto type = parse_type(m.value["t"].GetString());

        Id j=params->addProp(m.name.GetString(), type);

        auto const& valm = m.value.FindMember("v");
        if (valm == m.value.MemberEnd()) {
            switch (type) {
                case msys::IntType:
                    for (Id i=0; i<nparams; i++) {
                        params->value(i,j) = 0;
                    }
                    break;
                case msys::FloatType:
                    for (Id i=0; i<nparams; i++) {
                        params->value(i,j) = 0.0;
                    }
                    break;
                case msys::StringType:
                    for (Id i=0; i<nparams; i++) {
                        params->value(i,j) = "";
                    }
                    break;
            }
        } else {
            auto vals = m.value["v"].GetArray();
            CHECK_SIZE(vals, "v", nparams);
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
                        params->value(i,j) = get_name(names, vals[i].GetInt());
                    }
                    break;
            }
        }
    }
}

static void read_aux(Document const& d, SystemPtr mol) {
    auto const& tables = d.FindMember("aux");
    if (tables == d.MemberEnd()) return;
    auto& names = d["names"];
    for (auto& m : tables->value.GetObject()) {
        auto params = ParamTable::create();
        read_params(m.value, names, params, "aux params");
        mol->addAuxTable(m.name.GetString(), params);
    }
}

static void read_tables(Document const& d, SystemPtr mol) {
    auto const& names = d.FindMember("names");
    auto const& tables = d.FindMember("t");
    if (tables == d.MemberEnd()) return;
    for (auto& m : tables->value.GetObject()) {
        auto table = mol->addTable(m.name.GetString(), m.value["n"].GetInt());
        auto attrs = m.value.FindMember("a");
        if (attrs != m.value.MemberEnd()) {
            auto const& category = attrs->value.FindMember("c");
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
        read_params(m.value["p"], names->value, table->params(), table->name());
        auto& terms = m.value["t"];

        auto& particles = terms["i"];
        auto& params = terms["p"];
        CHECK_SIZE(params, "p", particles.Size() / table->atomCount());

        msys::IdList atoms(table->atomCount());
        Id particle_index = 0;
        for (Id i=0, n=params.Size(); i<n; i++) {
            for (Id j=0, m=table->atomCount(); j<m; j++) {
                atoms[j] = particles[particle_index++].GetInt();
            }
            table->addTerm(atoms, params[i].GetInt());
        }

        auto const &tags = m.value.FindMember("tags");
        if (tags != m.value.MemberEnd()) {
            read_tags(tags->value, table->props(), names->value);
        }
    }
}

static void read_cts(Document const& d, SystemPtr mol) {
    auto const& cts = d.FindMember("c");
    if (cts == d.MemberEnd()) return;

    Id ctid = 0;
    for (auto& ct : cts->value.GetArray()) {
        if (!ct.IsObject()) MSYS_FAIL("Object not found in ct array!");

        if (ctid != 0) {
            Id newctid = mol->addCt();
            if (ctid != newctid) MSYS_FAIL("Unexpected ct id encountered!");
        }
        auto &newct = mol->ct(ctid);

        // ct name
        auto const &name = ct.FindMember("n");
        if (name != ct.MemberEnd()) newct.setName(name->value.GetString());
           
        // ct key-values
        auto const &kv = ct.FindMember("k");
        if (kv != ct.MemberEnd()) {
            for (auto& m : kv->value.GetObject()) {
                Id id = newct.add(m.name.GetString(), msys::ValueType::StringType);
                newct.value(id) = m.value.GetString();
            }
        }

        // chain ct associations
        auto const &chains = ct.FindMember("c");
        if (chains != ct.MemberEnd()) {
            for (auto &chain_id : chains->value.GetArray()) {
                mol->setChain(chain_id.GetInt(), ctid);
            }
        }

        ctid += 1;
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
    read_cts(d, mol);
    msys::Analyze(mol);
    return mol;
}

namespace desres { namespace msys {

    SystemPtr ImportJson(std::string const& path) {
        Document d;
        auto json = slurp(path.data());
        ParseResult ok = d.Parse<kParseFullPrecisionFlag>(json.get());
        if (!ok) {
            MSYS_FAIL("Failed to parse JSON at " << ok.Offset() << ":" << GetParseError_En(ok.Code()));
            return nullptr;
        }
        SystemPtr mol = import_json(d);
        return mol;
    }

    SystemPtr ParseJson(const char* text) {
        Document d;
        ParseResult ok = d.Parse<kParseFullPrecisionFlag>(text);
        if (!ok) {
            MSYS_FAIL("Failed to parse JSON at " << ok.Offset() << ":" << GetParseError_En(ok.Code()));
            return nullptr;
        }
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
