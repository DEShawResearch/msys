#include "../json.hxx"

#include <rapidjson/document.h>
#include <rapidjson/writer.h>
#include <rapidjson/prettywriter.h>
#include <rapidjson/filewritestream.h>
#include <rapidjson/stringbuffer.h>

#include <unordered_map>
#include <fstream>
#include "../MsysThreeRoe.hpp"

using namespace rapidjson;
using namespace desres;
using msys::System;
using msys::Provenance;
using msys::SmallString;
using msys::TermTablePtr;
using msys::Id;

/* mapping from hash of string to index in document's names array */
typedef std::unordered_map<uint64_t, uint64_t> NameMap;
static uint64_t register_name(const char* ptr, size_t sz, NameMap& map, Document& d) {
    ThreeRoe tr;
    tr.Update(ptr, sz);
    auto id = tr.Final().first;
    auto pair = map.emplace(id, 0);
    if (pair.second) {
        auto& names = d["names"];
        auto& iter = pair.first;
        auto& alloc = d.GetAllocator();
        iter->second = names.Size();
        Value s;
        s.SetString(ptr, sz, alloc);
        names.PushBack(s, alloc);
    }
    return pair.first->second;
}

template <typename T>
static uint64_t register_name(T const& s, NameMap& map, Document& d) {
    return register_name(s.c_str(), s.size(), map, d);
}

static Value export_tags(msys::ParamTablePtr params, Document& d, NameMap& map) {
    auto& alloc = d.GetAllocator();
    Value tags(kObjectType);

    for (Id i=0, n=params->propCount(); i<n; i++) {
        Value tag(kObjectType);
        Value vals(kArrayType);
        Value ids(kArrayType);

        switch (params->propType(i)) {
            case msys::IntType:
                tag.AddMember("type", "i", alloc);
                for (Id j=0, m=params->paramCount(); j<m; j++) {
                    auto v = params->value(j,i).asInt();
                    if (v!=0) {
                        vals.PushBack(v, alloc);
                        ids.PushBack(j, alloc);
                    }
                }
                break;
            case msys::FloatType:
                tag.AddMember("type", "f", alloc);
                for (Id j=0, m=params->paramCount(); j<m; j++) {
                    auto v = params->value(j,i).asFloat();
                    if (v!=0) {
                        vals.PushBack(v, alloc);
                        ids.PushBack(j, alloc);
                    }
                }
                break;
            case msys::StringType:
                tag.AddMember("type", "s", alloc);
                for (Id j=0, m=params->paramCount(); j<m; j++) {
                    auto v = params->value(j,i).c_str();
                    if (*v) {
                        vals.PushBack(register_name(v, strlen(v), map, d), alloc);
                        ids.PushBack(j, alloc);
                    }
                }
                break;
        };
        if (true) { // vals.Size() > 0) {
            tag.AddMember("count", ids.Size(), alloc); 
            tag.AddMember("ids", ids, alloc);
            tag.AddMember("vals", vals, alloc);

            Value s;
            auto name = params->propName(i);
            s.SetString(name.data(), name.size(), alloc);
            tags.AddMember(s, tag, alloc);
        }
    }
    return tags;
}

static Value export_params(msys::ParamTablePtr params, Document& d, NameMap& map) {
    auto& alloc = d.GetAllocator();
    Value param(kObjectType);
    Value props(kObjectType);
    param.AddMember("count", params->paramCount(), alloc);

    for (Id i=0, n=params->propCount(); i<n; i++) {
        Value prop(kObjectType);
        Value vals(kArrayType);

        switch (params->propType(i)) {
            case msys::IntType:
                prop.AddMember("type", "i", alloc);
                for (Id j=0, m=params->paramCount(); j<m; j++) {
                    vals.PushBack(params->value(j,i).asInt(), alloc);
                }
                break;
            case msys::FloatType:
                prop.AddMember("type", "f", alloc);
                for (Id j=0, m=params->paramCount(); j<m; j++) {
                    vals.PushBack(params->value(j,i).asFloat(), alloc);
                }
                break;
            case msys::StringType:
                prop.AddMember("type", "s", alloc);
                for (Id j=0, m=params->paramCount(); j<m; j++) {
                    auto s = params->value(j,i).c_str();
                    vals.PushBack(register_name(s, strlen(s), map, d), alloc);
                }
                break;
        };
        prop.AddMember("vals", vals, alloc);

        Value s;
        auto name = params->propName(i);
        s.SetString(name.data(), name.size(), alloc);
        props.AddMember(s, prop, alloc);
    }
    param.AddMember("props", props, alloc);
    return param;
}

static Value export_cell(Document& d, NameMap& map, System& mol) {
    auto& alloc = d.GetAllocator();
    const double* cell = mol.global_cell[0];
    Value arr(kArrayType);
    for (int i=0; i<9; i++) {
        arr.PushBack(cell[i], alloc);
    }
    return arr;
}

static Value export_particles(Document& d, NameMap& map, System& mol) {
    auto& alloc = d.GetAllocator();
    Value p(kObjectType);
    Value names(kArrayType);
    Value anums(kArrayType);
    Value residues(kArrayType);
    Value x(kArrayType);
    Value y(kArrayType);
    Value z(kArrayType);
    Value vx(kArrayType);
    Value vy(kArrayType);
    Value vz(kArrayType);
    Value fc(kArrayType);
    Value mass(kArrayType);
    Value charge(kArrayType);

    for (auto i=mol.atomBegin(), e=mol.atomEnd(); i!=e; ++i) {
        auto const& a = mol.atomFAST(*i);
        names.PushBack(register_name(a.name,map,d), alloc);
        anums.PushBack(a.atomic_number, alloc);
        residues.PushBack(a.residue, alloc);
        x.PushBack(a.x, alloc);
        y.PushBack(a.y, alloc);
        z.PushBack(a.z, alloc);
        vx.PushBack(a.vx, alloc);
        vy.PushBack(a.vy, alloc);
        vz.PushBack(a.vz, alloc);
        fc.PushBack(a.formal_charge, alloc);
        mass.PushBack(a.mass, alloc);
        charge.PushBack(a.charge, alloc);
    }
    p.AddMember("count", mol.atomCount(), alloc);
    p.AddMember("name", names, alloc);
    p.AddMember("atomic_number", anums, alloc);
    p.AddMember("residue", residues, alloc);
    p.AddMember("formal_charge", fc, alloc);
    p.AddMember("mass", mass, alloc);
    p.AddMember("charge", charge, alloc);
    p.AddMember("x", x, alloc);
    p.AddMember("y", y, alloc);
    p.AddMember("z", z, alloc);
    p.AddMember("vx", vx, alloc);
    p.AddMember("vy", vy, alloc);
    p.AddMember("vz", vz, alloc);
    p.AddMember("tags", export_tags(mol.atomProps(), d, map), alloc);
    return p;
}

static Value export_bonds(Document& d, NameMap& map, System& mol) {
    auto& alloc = d.GetAllocator();
    Value b(kObjectType);
    Value p(kArrayType);
    Value o(kArrayType);
    for (auto i=mol.bondBegin(), e=mol.bondEnd(); i!=e; ++i) {
        auto const& a = mol.bondFAST(*i);
        p.PushBack(a.i, alloc);
        p.PushBack(a.j, alloc);
        o.PushBack(a.order, alloc);
    }
    b.AddMember("count", o.Size(), alloc);
    b.AddMember("particles", p, alloc);
    b.AddMember("order", o, alloc);
    b.AddMember("tags", export_tags(mol.bondProps(), d, map), alloc);
    return b;
}

static Value export_residues(Document& d, NameMap& map, System const& mol) {
    auto& alloc = d.GetAllocator();
    Value residues(kObjectType);
    Value chain(kArrayType);
    Value resid(kArrayType);
    Value name(kArrayType);
    Value insertion(kArrayType);
    for (auto i=mol.residueBegin(), e=mol.residueEnd(); i!=e; ++i) {
        auto const& a = mol.residueFAST(*i);
        chain.PushBack(a.chain, alloc);
        resid.PushBack(a.resid, alloc);
        name.PushBack(register_name(a.name,map,d), alloc);
        insertion.PushBack(register_name(a.insertion,map,d), alloc);
    }
    residues.AddMember("count", mol.residueCount(), alloc);
    residues.AddMember("chain", chain, alloc);
    residues.AddMember("resid", resid, alloc);
    residues.AddMember("name", name, alloc);
    // Make this optional
    residues.AddMember("insertion", insertion, alloc);
    return residues;
}

static Value export_chains(Document& d, NameMap& map, System const& mol) {
    auto& alloc = d.GetAllocator();
    Value chains(kObjectType);
    Value name(kArrayType);
    Value segid(kArrayType);
    for (auto i=mol.chainBegin(), e=mol.chainEnd(); i!=e; ++i) {
        auto const& a = mol.chainFAST(*i);
        name.PushBack(register_name(a.name,map,d), alloc);
        segid.PushBack(register_name(a.segid,map,d), alloc);
    }
    chains.AddMember("name", name, alloc);
    // Make this optional
    chains.AddMember("segid", segid, alloc);
    chains.AddMember("count", mol.chainCount(), alloc);
    return chains;
}

static void export_terms(Value& obj, TermTablePtr table, Document& d) {
    auto& alloc = d.GetAllocator();
    Value terms(kObjectType);
    Value particles(kArrayType);
    Value params(kArrayType);

    for (auto const& term : *table) {
        for (Id i=0, n=table->atomCount(); i<n; i++) {
            particles.PushBack(term.atom(i), alloc);
        }
        params.PushBack(int32_t(term.param()), alloc);
    }
    terms.AddMember("count", table->termCount(), alloc);
    terms.AddMember("particles", particles, alloc);
    terms.AddMember("params", params, alloc);
    obj.AddMember("arity", table->atomCount(), alloc);
    obj.AddMember("terms", terms, alloc);
}

static Value export_tables(Document& d, NameMap& map, System const& mol) {
    auto& alloc = d.GetAllocator();
    Value tables(kObjectType);
    for (auto table_name : mol.tableNames()) {
        auto table = mol.table(table_name);
        Value tableobj(kObjectType);
        Value s;

        export_terms(tableobj, table, d);
        tableobj.AddMember("param", export_params(table->params(), d, map), alloc);
        tableobj.AddMember("tags", export_tags(table->props(), d, map), alloc);

        Value attrs(kObjectType);
        auto category = msys::print(table->category);
        s.SetString(category.data(), category.size(), alloc);
        attrs.AddMember("category", s, alloc);
        if (table_name == "nonbonded") {
            auto rule = mol.nonbonded_info.vdw_rule;
            s.SetString(rule.data(), rule.size(), alloc);
            attrs.AddMember("vdw_rule", s, alloc);
        }
        tableobj.AddMember("attrs", attrs, alloc);

        s.SetString(table_name.data(), table_name.size(), alloc);
        tables.AddMember(s, tableobj, alloc);
    }
    return tables;
}

static void export_json(Document& d, System& mol, Provenance const& provenance, unsigned flags) {

    auto& alloc = d.GetAllocator();
    d.SetObject();
    Value names(kArrayType);
    d.AddMember("names", names, alloc);

    NameMap map;
    d.AddMember("cell", export_cell(d, map, mol), alloc);
    d.AddMember("particles", export_particles(d, map, mol), alloc);
    d.AddMember("bonds", export_bonds(d, map, mol), alloc);
    d.AddMember("residues", export_residues(d, map, mol), alloc);
    d.AddMember("chains", export_chains(d, map, mol), alloc);
    if (!(flags & desres::msys::JsonExport::StructureOnly)) {
        d.AddMember("tables", export_tables(d, map, mol), alloc);
    }
}

namespace desres { namespace msys {

    void ExportJson(SystemPtr mol, std::string const& path, Provenance const& provenance,
            unsigned flags) {

        char writeBuffer[65536];
        auto fp = std::shared_ptr<FILE>(fopen(path.data(), "wb"), fclose);
        if (!fp) {
            MSYS_FAIL(strerror(errno));
        }

        Document d;
        export_json(d, *mol, provenance, flags);

        FileWriteStream os(fp.get(), writeBuffer, sizeof(writeBuffer));
        if (flags & JsonExport::Whitespace) {
            PrettyWriter<FileWriteStream> writer(os);
            d.Accept(writer);
        } else {
            Writer<FileWriteStream> writer(os);
            d.Accept(writer);
        }
    }

    std::string FormatJson(SystemPtr mol, Provenance const& provenance, unsigned flags) {

        Document d;
        export_json(d, *mol, provenance, flags);
        StringBuffer buffer;
        if (flags & JsonExport::Whitespace) {
            PrettyWriter<StringBuffer> writer(buffer);
            d.Accept(writer);
        } else {
            Writer<StringBuffer> writer(buffer);
            d.Accept(writer);
        }
        return buffer.GetString();
    }

}}

