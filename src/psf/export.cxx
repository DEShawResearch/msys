#include "../psf.hxx"
#include "../clone.hxx"

using namespace desres::msys;

void desres::msys::ExportPSF(SystemPtr mol, std::string const& path) {
    FILE* fp = fopen(path.data(), "w");
    if (!fp) {
        MSYS_FAIL("Unable to open '" << path << "' for writing");
    }
    std::shared_ptr<FILE> dtor(fp, fclose);

    if (mol->atomCount() != mol->maxAtomId() ||
        mol->bondCount() != mol->maxBondId()) {
        mol = Clone(mol, mol->atoms());
    }

    fprintf(fp, "PSF\n\n%8d !NTITLE\n", 1);
    fprintf(fp, "REMARKS %s\n", "msys generated charmm psf file");
    fprintf(fp, "\n");

    fprintf(fp, "%8d !NATOM\n", mol->atomCount());
    for (auto id : mol->atoms()) {
        auto& atm = mol->atomFAST(id);
        auto& res = mol->residueFAST(atm.residue);
        auto& chn = mol->chainFAST(res.chain);

        fprintf(fp, "%8d %-4s %-4d %-4s %-4s %4d %10.6f     %9.4f  %10d\n",
                id+1, chn.segid.c_str(), res.resid, res.name.c_str(),
                atm.name.c_str(), 0, atm.charge, atm.mass, 0);
    }
    fprintf(fp, "\n");

    fprintf(fp, "%8d !NBOND: bonds\n", mol->bondCount());
    for (auto id : mol->bonds()) {
        auto& bnd = mol->bondFAST(id);
        fprintf(fp, "%8d%8d", bnd.i+1, bnd.j+1);
        if ((id % 4) == 3) {
            fprintf(fp, "\n");
        }
    }
    if (mol->bondCount() % 4 != 0) {
        fprintf(fp, "\n");
    }

    fprintf(fp, "\n");
}


