#include "view.hxx"
#include "import_dms.hxx"

using namespace desres::msys;

int main(int argc, char *argv[]) {
    System msys;
    for (int i=1; i<argc; i++) {
        desres::msys::ImportDMS(msys, argv[i]);
        Id n=msys.atomCount();
        IdList atoms(n/10);
        for (Id i=0; i<n/10; i++) atoms[i]=i;
        View view(&msys, atoms);
        assert(view.atoms().size()==n/10);
    }
    return 0;
}
