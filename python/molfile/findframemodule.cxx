#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <msys/molfile/findframe.hxx>
#include <msys/types.hxx>

namespace py = pybind11;
namespace ff = desres::molfile::findframe;

namespace {
    struct Oracle {
        const double * times;
        Oracle(const double *_times) : times(_times) {}
        double operator[](ssize_t i) const { return times[i]; }
    };

    typedef ssize_t (*findfunc)(ssize_t N, const Oracle& times, double T);

    template <findfunc f> 
    ssize_t wrap(py::array_t<double, py::array::c_style | py::array::forcecast> arr, double  T) {
        auto times = arr.data();
        auto N = arr.size();
        return f(N, Oracle(times), T);
    }
}

PYBIND11_MODULE(findframe, m) {
    m.def("at_time_near", wrap<ff::at_time_near<Oracle> >);
    m.def("at_time_lt", wrap<ff::at_time_lt<Oracle> >);
    m.def("at_time_le", wrap<ff::at_time_le<Oracle> >);
    m.def("at_time_gt", wrap<ff::at_time_gt<Oracle> >);
    m.def("at_time_ge", wrap<ff::at_time_ge<Oracle> >);
}

#if _GLIBCXX_RELEASE >= 11
namespace std {

/* This avoids the GLIBCXX_3.4.29 symbol version. */
void __attribute__((weak)) __throw_bad_array_new_length() { MSYS_FAIL("bad array new length"); }

}  // namespace std
#endif
