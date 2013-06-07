#include "atommatch.hxx"
#include <boost/python.hpp>
#include <msys/version.hxx>

using namespace desres;
using namespace desres::fep_atommatch;
using namespace boost::python;

namespace {

    class ScoreFctPy : public ScoreFct {
        public:
            static ScoreFctPtr create(const object& py_apply) {
                return ScoreFctPtr(new ScoreFctPy(py_apply));
            }
            virtual double apply(msys::SystemPtr mol1, const msys::IdList& atoms1,
                    msys::SystemPtr mol2, const msys::IdList& atoms2) const {
                list L1;
                list L2;
                BOOST_FOREACH(msys::Id id, atoms1) { L1.append(id); }
                BOOST_FOREACH(msys::Id id, atoms2) { L2.append(id); }
                double output;
                try {
                    output = extract<double>(_py_apply(mol1, L1, mol2, L2));
                } catch (std::exception& e) {
                    FAIL("Error calling Python function: " + std::string(e.what()));
                }
                return output;
            }
        private:
            ScoreFctPy(const object& py_apply) : _py_apply(py_apply) { }
            object _py_apply;
    };

    list fep_atom_match(msys::SystemPtr mol1, msys::SystemPtr mol2, ScoreFctPtr rep) {
        MatchList match;
        FEPAtomMatch(mol1, mol2, rep, match);
        list L;
        for (unsigned i = 0; i < match.size(); ++i)
            L.append(make_tuple(match[i].first, match[i].second));
        return L;
    }

}

BOOST_PYTHON_MODULE(_atommatch) {
    boost::python::def("FEPAtomMatch", fep_atom_match);
    class_<ScoreFct, ScoreFctPtr, boost::noncopyable>("ScoreFct", no_init)
        .def("__init__", make_constructor(&ScoreFctPy::create))
        .def("apply", &ScoreFctPy::apply)
        ;
}
