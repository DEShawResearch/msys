#include "atommatch.hxx"
#include <boost/python.hpp>
#include <msys/version.hxx>

using namespace desres::msys;
using namespace desres::msys::atommatch;
using namespace boost::python;

namespace {

    class ScoreFctPy : public ScoreFct {
        public:
            static ScoreFctPtr create(const object& py_apply) {
                return ScoreFctPtr(new ScoreFctPy(py_apply));
            }
            virtual double apply(SystemPtr mol1, const IdList& atoms1,
                    SystemPtr mol2, const IdList& atoms2) const {
                list L1;
                list L2;
                BOOST_FOREACH(Id id, atoms1) { L1.append(id); }
                BOOST_FOREACH(Id id, atoms2) { L2.append(id); }
                double output;
                try {
                    output = extract<double>(_py_apply(mol1, L1, mol2, L2));
                } catch (std::exception& e) {
                    MSYS_FAIL("Error calling Python function: " + std::string(e.what()));
                }
                return output;
            }
        private:
            ScoreFctPy(const object& py_apply) : _py_apply(py_apply) { }
            object _py_apply;
    };

    list atom_match(SystemPtr mol1, SystemPtr mol2, ScoreFctPtr rep) {
        MatchList match;
        AtomMatch(mol1, mol2, rep, match);
        list L;
        for (unsigned i = 0; i < match.size(); ++i)
            L.append(make_tuple(match[i].first, match[i].second));
        return L;
    }

}

BOOST_PYTHON_MODULE(_atommatch) {
    boost::python::def("AtomMatch", atom_match);
    class_<ScoreFct, ScoreFctPtr, boost::noncopyable>("ScoreFct", no_init)
        .def("__init__", make_constructor(&ScoreFctPy::create))
        .def("apply", &ScoreFctPy::apply)
        ;
}
