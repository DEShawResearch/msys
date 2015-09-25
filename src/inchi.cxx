#include "inchi.hxx"
#include "clone.hxx"
#include "elements.hxx"
#include "analyze.hxx"
#include <string.h>

using namespace desres::msys;

#ifdef MSYS_WITHOUT_INCHI
InChI InChI::create(SystemPtr mol, unsigned options) {
    MSYS_FAIL("No InChI support");
}
String InChI::key() const {
    MSYS_FAIL("No InChI support");
}
#else

#include <inchi/inchi_api.h>

std::vector<InChI> InChI::analyze(SystemPtr mol, unsigned options) {
    MultiIdList fragments;
    mol->updateFragids(&fragments);
    IdList distinct = FindDistinctFragments(mol, fragments);
    std::vector<InChI> result;
    for (Id id : distinct) {
        if (fragments[id].size()<1024) {
            result.emplace_back(create(Clone(mol, fragments[id]), options));
        }
    }
    return result;
}

InChI InChI::create(SystemPtr mol, unsigned options) {

    if (mol->atomCount() >= 1024) {
        MSYS_FAIL("Too many atoms: " << mol->atomCount() << " >= 1024");
    }

    /* make the atoms contiguous */
    if (mol->maxAtomId() != mol->atomCount()) {
        mol = Clone(mol, mol->atoms());
    }

    inchi_Input  input[1];
    inchi_Output output[1];
    memset(input,0,sizeof(input));
    memset(output,0,sizeof(output));

    std::vector<inchi_Atom> atoms(mol->atomCount());

    for (Id i=0; i<mol->maxAtomId(); i++) {
        inchi_Atom& iatom = atoms[i];
        atom_t const& matom = mol->atom(i);

        iatom.x = matom.x;
        iatom.y = matom.y;
        iatom.z = matom.z;
        iatom.num_bonds = mol->bondCountForAtom(i);
        Id j=0;
        for (Id b : mol->bondsForAtom(i)) {
            bond_t const& bond = mol->bond(b);
            iatom.neighbor[j] = bond.other(i);
            iatom.bond_type[j]= bond.order;
            ++j;
        }
        strcpy(iatom.elname, AbbreviationForElement(matom.atomic_number));
        iatom.charge = matom.formal_charge;
    }

    if (atoms.size()) input->atom = &atoms[0];
    input->num_atoms = atoms.size();

    std::string ioptions;
    if (options) {
        if (options & DoNotAddH) ioptions += " -DoNotAddH";
        if (options & SNon)      ioptions += " -SNon";
        if (options & FixedH)    ioptions += " -FixedH";

        input->szOptions = &ioptions[0];
    }

    int rc = GetINCHI(input, output);
    boost::shared_ptr<inchi_Output> dtor(output, FreeINCHI);

    if (!(rc==0 || rc==1)) {
        MSYS_FAIL("InChI failed: " << output->szMessage);
    }
    return InChI(rc, output->szInChI, output->szAuxInfo, output->szMessage);
}

String InChI::key() const {
    char buf[32];
    GetINCHIKeyFromINCHI(_string.c_str(), 0,0, buf, NULL, NULL);
    return buf;
}
#endif

