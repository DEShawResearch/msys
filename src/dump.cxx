#include "dump.hxx"
#include <fastjson/print.hxx>
#include <stdio.h>

using namespace desres::msys;
using desres::fastjson::floatify;

namespace {
    class Dump {
        FILE* fd;
        const unsigned flags;
        char buf[18];

        const char* floated(double v) {
            floatify(v,buf);
            return buf;
        }

    public:
        Dump(FILE* fd, unsigned flags) 
        : fd(fd), flags(flags) {}
        void highlevel(SystemPtr mol);
        void cell(SystemPtr mol);
        void atoms(SystemPtr mol);
        void bonds(SystemPtr mol);
        void residues(SystemPtr mol);
        void chains(SystemPtr mol);
        void table(TermTablePtr );
        void tables(SystemPtr mol);
        
        void dump(SystemPtr mol);
    };
}

void Dump::highlevel(SystemPtr mol) {
    fprintf(fd, "System:   %s\n", mol->name.c_str());
    fprintf(fd, "Atoms:    %u\n", mol->atomCount());
    fprintf(fd, "Bonds:    %u\n", mol->bondCount());
    fprintf(fd, "Residues: %u\n", mol->residueCount());
    fprintf(fd, "Chains:   %u\n", mol->chainCount());
    fprintf(fd, "Cts:      %u\n", mol->ctCount());
}

void Dump::cell(SystemPtr mol) {
    GlobalCell const& cell = mol->global_cell;
    for (int i=0; i<3; i++) {
        fprintf(fd, "Cell %c:   ", 'A'+i);
        for (int j=0; j<3; j++) fprintf(fd, "%s ", floated(cell[i][j]));
        fprintf(fd, "\n");
    }
}

void Dump::atoms(SystemPtr mol) {
    IdList const& atoms = mol->atoms();
    for (Id i=0; i<atoms.size(); i++) {
        Id atm = atoms[i];
        atom_t const& atom = mol->atom(atm);
        fprintf(fd, "Atom %-5u fragid %-5u residue %-5u bonds %-2u anum %-2d fq %-2d name %s\n",
                atm, atom.fragid, atom.residue, mol->bondCountForAtom(atm),
                atom.atomic_number, atom.formal_charge, atom.name.c_str());
    }
}

void Dump::bonds(SystemPtr mol) {
    IdList const& bonds = mol->bonds();
    for (Id i=0; i<bonds.size(); i++) {
        Id bnd = bonds[i];
        bond_t const& bond = mol->bond(bnd);
        fprintf(fd, "Bond %-5u ai %-5u aj %-5u order %d\n", 
                bnd, bond.i, bond.j, bond.order);
    }
}

void Dump::residues(SystemPtr mol) {
    IdList const& list = mol->residues();
    for (Id i=0; i<list.size(); i++) {
        Id id = list[i];
        residue_t const& res = mol->residue(id);
        fprintf(fd, "Residue %-5u chain %-5u atoms %-4u resid %-6d name %s",
            id, res.chain, mol->atomCountForResidue(id), res.resid,
            res.name.c_str());
        if (res.insertion.size()) fprintf(fd, " (%s)", res.insertion.c_str());
        fprintf(fd, "\n");
    }
}

void Dump::chains(SystemPtr mol) {
    IdList const& list = mol->chains();
    for (Id i=0; i<list.size(); i++) {
        Id id = list[i];
        chain_t const& chn = mol->chain(id);
        fprintf(fd, "Chain %-5u ct %-6u residues %-5u name %-4s",
                id, chn.ct, mol->residueCountForChain(id), chn.name.c_str());
        if (chn.segid.size()) fprintf(fd, " segid %s", chn.segid.c_str());
        fprintf(fd, "\n");
    }
}

void Dump::table(TermTablePtr table) {
    ParamTablePtr params = table->params();
    const Id natoms = table->atomCount();
    fprintf(fd, "Table %-26s atoms %-2u terms %-6u params %-6u props %u\n",
            table->name().c_str(), table->atomCount(), table->termCount(),
            params->paramCount(), params->propCount());
    std::vector<String> termprops, paramprops;
    IdList termpropids, parampropids, termproptypes, paramproptypes;
    std::set<std::string> paraminfo;
    paraminfo.insert("ff");
    paraminfo.insert("memo");
    paraminfo.insert("type");
    paraminfo.insert("comment");
    paraminfo.insert("typekey");
    for (Id i=0; i<table->termPropCount(); i++) {
        termprops.push_back(table->termPropName(i));
    }
    for (Id i=0; i<params->propCount(); i++) {
        String const& name = params->propName(i);
        if (flags & TextExport::ParamInfo || !paraminfo.count(name)) {
            paramprops.push_back(name);
        }
    }
    std::sort(termprops.begin(), termprops.end());
    std::sort(paramprops.begin(), paramprops.end());
    for (Id i=0; i<termprops.size(); i++) {
        termpropids.push_back(table->termPropIndex(termprops[i]));
        termproptypes.push_back(table->termPropType(termpropids[i]));
    }
    for (Id i=0; i<paramprops.size(); i++) {
        parampropids.push_back(params->propIndex(paramprops[i]));
        paramproptypes.push_back(params->propType(parampropids[i]));
    }
    IdList const& terms = table->terms();

    for (Id i=0; i<terms.size(); i++) {
        Id term = terms[i];
        fprintf(fd, "Term %-6u atoms ", term);
        for (Id j=0; j<natoms; j++) fprintf(fd,"%-6u ", table->atom(term,j));
        for (Id j=0; j<paramprops.size(); j++) {
            ValueRef v = table->propValue(term,parampropids[j]);
            std::string const& nm = paramprops[j];
            const char* n = nm.c_str();
            switch (paramproptypes[j]) {
                case IntType:   fprintf(fd,"%s %-8ld ", n, v.asInt()); break;
                case FloatType: fprintf(fd,"%s %-8s ", n, floated(v.asFloat())); break;
                default:        fprintf(fd,"%s %-8s ", n, v.asString().c_str()); break;
            }
        }
        for (Id j=0; j<termprops.size(); j++) {
            ValueRef v = table->termPropValue(term, termpropids[j]);
            std::string const& nm = termprops[j];
            const char* n = nm.c_str();
            switch (termproptypes[j]) {
                case IntType:   fprintf(fd,"%s %-8ld ", n, v.asInt()); break;
                case FloatType: fprintf(fd,"%s %-8s ", n, floated(v.asFloat())); break;
                default:        fprintf(fd,"%s %-8s ", n, v.asString().c_str()); break;
            }
        }
        fprintf(fd,"\n");
    }

}

void Dump::tables(SystemPtr mol) {
    std::vector<String> names = mol->tableNames();
    for (Id i=0; i<names.size(); i++) {
        table(mol->table(names[i]));
    }
}

void Dump::dump(SystemPtr mol) {
    highlevel(mol);
    cell(mol);
    atoms(mol);
    bonds(mol);
    residues(mol);
    chains(mol);
    tables(mol);
}

void desres::msys::ExportText(SystemPtr mol, FILE* fp,
                              unsigned flags) {
    Dump(fp,flags).dump(mol);
}
