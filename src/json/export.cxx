#include "../json.hxx"

#if defined __has_include
#  if __has_include (<rapidjson/document.h>)
#include <rapidjson/document.h>
#include <rapidjson/writer.h>
#include <rapidjson/prettywriter.h>
#include <rapidjson/filewritestream.h>
#include <rapidjson/stringbuffer.h>
#    define MSYS_WITH_RAPID_JSON
using namespace rapidjson;
#  endif
#endif




#include <unordered_map>
#include <fstream>
#include "../MsysThreeRoe.hpp"

using namespace desres;
using msys::System;
using msys::Provenance;
using msys::SmallString;
using msys::TermTablePtr;
using msys::Id;

#if defined MSYS_WITH_RAPID_JSON

/* mapping from hash of string to index in document's names array */
typedef std::unordered_map<uint64_t, uint64_t> NameMap;
static uint64_t register_name(const char* ptr, size_t sz, NameMap& map, Document& d) {
    // Elide empty strings from the JSON output
    if (sz == 0)
        return 0;

    // Then need to register the empty string as the first name if any names are given
    if (d.FindMember("names") == d.MemberEnd()) {
        Value names(kArrayType);
        auto& alloc = d.GetAllocator();
        Value s;
        s.SetString("", 0, alloc);
        names.PushBack(s, alloc);
        d.AddMember("names", names, alloc);
    }

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
    bool nonzero = false;
    for (int i=0; i<9; i++) {
        arr.PushBack(cell[i], alloc);
        if (cell[i] != 0) nonzero = true;
    }
    if (!nonzero) arr.Clear();
    return arr;
}

static Value export_particles(Document& d, NameMap& map, System& mol) {
    auto& alloc = d.GetAllocator();
    Value p(kObjectType);
    Value names(kArrayType);
    Value anums(kArrayType);
    Value residues(kArrayType);
    Value pos(kArrayType);
    Value vel(kArrayType);
    Value fc(kArrayType);
    Value mass(kArrayType);
    Value charge(kArrayType);

    bool nonzero_nam=false;
    bool nonzero_anm=false;
    bool nonzero_fch=false;
    bool nonzero_pos=false;
    bool nonzero_vel=false;
    bool nonzero_res=false;

    for (auto i=mol.atomBegin(), e=mol.atomEnd(); i!=e; ++i) {
        auto const& a = mol.atomFAST(*i);
        uint64_t name = register_name(a.name, map, d);
        names.PushBack(name, alloc);
        anums.PushBack(a.atomic_number, alloc);
        residues.PushBack(a.residue, alloc);
        pos.PushBack(a.x, alloc);
        pos.PushBack(a.y, alloc);
        pos.PushBack(a.z, alloc);
        vel.PushBack(a.vx, alloc);
        vel.PushBack(a.vy, alloc);
        vel.PushBack(a.vz, alloc);
        fc.PushBack(a.formal_charge, alloc);
        mass.PushBack(a.mass, alloc);
        charge.PushBack(a.charge, alloc);
        if (name!=0) nonzero_nam=true;
        if (a.atomic_number!=0) nonzero_anm=true;
        if (a.formal_charge!=0) nonzero_fch=true;
        if (a.residue!=0) nonzero_res=true;
        if (a.x!=0  || a.y!=0  || a.z!=0)  nonzero_pos=true;
        if (a.vx!=0 || a.vy!=0 || a.vz!=0) nonzero_vel=true;
    }
    p.AddMember("count", mol.atomCount(), alloc);
    if (nonzero_nam) p.AddMember("name", names, alloc);
    if (nonzero_anm) p.AddMember("atomic_number", anums, alloc);
    if (nonzero_fch) p.AddMember("formal_charge", fc, alloc);
    p.AddMember("mass", mass, alloc);
    p.AddMember("charge", charge, alloc);
    if (nonzero_pos) p.AddMember("position", pos, alloc);
    if (nonzero_vel) p.AddMember("velocity", vel, alloc);
    if (nonzero_res) p.AddMember("residue", residues, alloc);
    
    Value tags = export_tags(mol.atomProps(), d, map);
    if (!tags.ObjectEmpty()) p.AddMember("tags", tags, alloc);
    return p;
}

static Value export_bonds(Document& d, NameMap& map, System& mol) {
    auto& alloc = d.GetAllocator();
    Value b(kObjectType);
    Value p(kArrayType);
    Value o(kArrayType);

    bool nonzero_order=false;

    for (auto i=mol.bondBegin(), e=mol.bondEnd(); i!=e; ++i) {
        auto const& a = mol.bondFAST(*i);
        p.PushBack(a.i, alloc);
        p.PushBack(a.j, alloc);
        o.PushBack(a.order, alloc);
        if (a.order!=0) nonzero_order=true;
    }
    b.AddMember("count", o.Size(), alloc);
    b.AddMember("particles", p, alloc);
    if (nonzero_order) b.AddMember("order", o, alloc);
    Value tags = export_tags(mol.bondProps(), d, map);
    if (!tags.ObjectEmpty()) b.AddMember("tags", tags, alloc);
    return b;
}

static Value export_residues(Document& d, NameMap& map, System const& mol) {
    auto& alloc = d.GetAllocator();
    Value residues(kObjectType);
    Value chain(kArrayType);
    Value resid(kArrayType);
    Value name(kArrayType);
    Value insertion(kArrayType);
    bool nonzero = false;
    for (auto i=mol.residueBegin(), e=mol.residueEnd(); i!=e; ++i) {
        auto const& a = mol.residueFAST(*i);
        chain.PushBack(a.chain, alloc);
        resid.PushBack(a.resid, alloc);
        name.PushBack(register_name(a.name,map,d), alloc);
        insertion.PushBack(register_name(a.insertion,map,d), alloc);
        if (a.chain!=0 || a.resid!=0 || !a.name.empty() || !a.insertion.empty()) {
            nonzero = true;
        }
    }
    if (nonzero || mol.residueCount() > 1) {
        residues.AddMember("count", mol.residueCount(), alloc);
        residues.AddMember("chain", chain, alloc);
        residues.AddMember("resid", resid, alloc);
        residues.AddMember("name", name, alloc);
        residues.AddMember("insertion", insertion, alloc);
    }
    return residues;
}

static Value export_chains(Document& d, NameMap& map, System const& mol) {
    auto& alloc = d.GetAllocator();
    Value chains(kObjectType);
    Value name(kArrayType);
    Value segid(kArrayType);
    bool nonzero = false;
    for (auto i=mol.chainBegin(), e=mol.chainEnd(); i!=e; ++i) {
        auto const& a = mol.chainFAST(*i);
        name.PushBack(register_name(a.name,map,d), alloc);
        segid.PushBack(register_name(a.segid,map,d), alloc);
        if (!(a.name.empty() && a.segid.empty())) nonzero=false;
    }
    if (nonzero || mol.chainCount()>1) {
        chains.AddMember("name", name, alloc);
        chains.AddMember("segid", segid, alloc);
        chains.AddMember("count", mol.chainCount(), alloc);
    }
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

static Value export_aux(Document& d, NameMap& map, System const& mol) {
    auto& alloc = d.GetAllocator();
    Value aux(kObjectType);
    for (auto& name : mol.auxTableNames()) {
        Value s;
        s.SetString(name.data(), name.size(), alloc);
        aux.AddMember(s, export_params(mol.auxTable(name), d, map), alloc);
    }
    return aux;
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
        Value tags = export_tags(table->props(), d, map);
        if (!tags.ObjectEmpty()) tableobj.AddMember("tags", tags, alloc);

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

    NameMap map;
    Value cell = export_cell(d, map, mol);
    if (cell.Size() != 0) d.AddMember("cell", cell, alloc);
    d.AddMember("particles", export_particles(d, map, mol), alloc);
    d.AddMember("bonds", export_bonds(d, map, mol), alloc);
    Value residues = export_residues(d, map, mol);
    if (!residues.ObjectEmpty()) d.AddMember("residues", residues, alloc);
    Value chains = export_chains(d, map, mol);
    if (!chains.ObjectEmpty()) d.AddMember("chains", chains, alloc);
    if (!(flags & desres::msys::JsonExport::StructureOnly)) {
        d.AddMember("tables", export_tables(d, map, mol), alloc);

        Value aux = export_aux(d, map, mol);
        if (!aux.ObjectEmpty()) d.AddMember("aux", aux, alloc);
    }
}

namespace desres { namespace msys {

    void ExportJson(SystemPtr mol, std::string const& path, Provenance const& provenance,
            unsigned flags, int maxDecimals) {

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
            if (maxDecimals >= 0) {
                writer.SetMaxDecimalPlaces(maxDecimals);
            }
            d.Accept(writer);
        } else {
            Writer<FileWriteStream> writer(os);
            if (maxDecimals >= 0) {
                writer.SetMaxDecimalPlaces(maxDecimals);
            }
            d.Accept(writer);
        }
    }

    std::string FormatJson(SystemPtr mol, Provenance const& provenance, unsigned flags, int maxDecimals) {

        Document d;
        export_json(d, *mol, provenance, flags);
        StringBuffer buffer;
        if (flags & JsonExport::Whitespace) {
            PrettyWriter<StringBuffer> writer(buffer);
            if (maxDecimals >= 0) {
                writer.SetMaxDecimalPlaces(maxDecimals);
            }
            d.Accept(writer);
        } else {
            Writer<StringBuffer> writer(buffer);
            if (maxDecimals >= 0) {
                writer.SetMaxDecimalPlaces(maxDecimals);
            }
            d.Accept(writer);
        }
        return buffer.GetString();
    }

}}

#else
#warning "rapidjson not available; no json export support"
namespace desres { namespace msys {

    void ExportJson(SystemPtr mol, std::string const& path, Provenance const& provenance,
            unsigned flags, int maxDecimals) {
        MSYS_FAIL("rapidjson not available; no json export support");
            }



    std::string FormatJson(SystemPtr mol, Provenance const& provenance, unsigned flags, int maxDecimals) {
         MSYS_FAIL("rapidjson not available; no json export support");
    }

}}

#endif
