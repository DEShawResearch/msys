#include <boost/python.hpp>
#include <numpy/arrayobject.h>

#include "findframe.hxx"

namespace bp = boost::python;
namespace ff = desres::molfile::findframe;

namespace {
    struct Oracle {
        const double * times;
        Oracle(const double *_times) : times(_times) {}
        double operator[](ssize_t i) const { return times[i]; }
    };

    typedef ssize_t (*findfunc)(ssize_t N, const Oracle& times, double T);

    template <findfunc f> 
    ssize_t wrap(bp::object& obj, double T) {
        PyObject * arr = PyArray_FROMANY(
                obj.ptr(),      /* PyObject */
                NPY_DOUBLE,     /* type */
                1, 1,           /* min, max depth (insist on 1d sequence) */
                NPY_IN_ARRAY);  /* get back contiguous and aligned */

        if (!arr) return -1;

        const double * times = (double *)PyArray_DATA(arr);
        ssize_t N = PyArray_DIM(arr,0);
        ssize_t result = f(N, Oracle(times), T);
        Py_DECREF(arr);
        return result;
    }
}

BOOST_PYTHON_MODULE(findframe) {
    import_array();

    bp::def("at_time_near", wrap<ff::at_time_near<Oracle> >);
    bp::def("at_time_lt", wrap<ff::at_time_lt<Oracle> >);
    bp::def("at_time_le", wrap<ff::at_time_le<Oracle> >);
    bp::def("at_time_gt", wrap<ff::at_time_gt<Oracle> >);
    bp::def("at_time_ge", wrap<ff::at_time_ge<Oracle> >);
}
