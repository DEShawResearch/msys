#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include "spatial_hash.hxx"

using namespace pybind11;
using namespace desres::msys;

using pos_t = array_t<float, array::forcecast | array::c_style>;
using box_t = array_t<double, array::forcecast | array::c_style>;
using ids_t = array_t<uint32_t, array::forcecast | array::c_style>;

template <typename Func, typename Param>
static object hash_find(SpatialHash& hash,
                           Param r,
                           pos_t posarr,
                           ids_t idsarr,
                           Func func) {

    if (posarr.ndim() != 2 || posarr.shape(1) != 3) {
        throw std::invalid_argument("expected Nx3 array for pos");
    }
    auto pos = posarr.data();
    auto ids = idsarr.data();
    auto n = idsarr.shape(0);
    uint32_t npos = posarr.shape(0);
    for (int i=0; i<n; i++) {
        if (ids[i] >= npos) {
            PyErr_Format(PyExc_ValueError, "index out of bounds: %d", ids[i]);
            throw error_already_set();
        }
    }
    IdList result = (hash.*func)(r, pos, n, ids);
    ssize_t sz = result.size();
    auto arr = ids_t(sz);
    if (!result.empty()) {
        memcpy(arr.mutable_data(), result.data(), result.size()*sizeof(result[0]));
    }
    return arr;
}

static object hash_find_within(SpatialHash& hash,
                                  float r,
                                  pos_t pos,
                                  ids_t ids,
                                  bool reuse_voxels) {

    return reuse_voxels ? hash_find(hash, r, pos, ids, &SpatialHash::find_within)
                        : hash_find(hash, r, pos, ids, &SpatialHash::findWithin);
}

static object hash_find_nearest(SpatialHash& hash,
                                   int k,
                                   pos_t pos,
                                   ids_t ids) {

    return hash_find(hash, k, pos, ids, &SpatialHash::findNearest);
}

static object hash_find_contacts(SpatialHash& hash,
                                    float r,
                                    pos_t posarr,
                                    object idsobj,
                                    bool reuse_voxels) {

    if (posarr.ndim() != 2 || posarr.shape(1) != 3) {
        throw std::invalid_argument("expected Nx3 array for pos");
    }
    auto pos = posarr.data();
    uint32_t npos = posarr.shape(0);

    const Id* ids = NULL;
    auto idsarr = ids_t::ensure(idsobj);
    if (idsarr) {
        ids = idsarr.data();
        auto n = idsarr.shape(0);
        for (int i=0; i<n; i++) {
            if (ids[i] >= npos) {
                PyErr_Format(PyExc_ValueError, "index out of bounds: %d", ids[i]);
                throw error_already_set();
            }
        }
        npos = n;
    }

    SpatialHash::contact_array_t contacts;

    if (reuse_voxels) {
        // drop the GIL here?
        hash.findContactsReuseVoxels(r, pos, npos, ids, &contacts);
    } else {
        hash.findContacts(r, pos, npos, ids, &contacts);
    }

    auto dim = contacts.count;
    std::for_each(contacts.d2, contacts.d2+dim, [](float& x) {x=std::sqrt(x);});

    auto iarr = ids_t(dim);
    auto jarr = ids_t(dim);
    auto darr = pos_t(dim);
    memcpy(iarr.mutable_data(), contacts.i, dim*sizeof(*contacts.i));
    memcpy(jarr.mutable_data(), contacts.j, dim*sizeof(*contacts.j));
    memcpy(darr.mutable_data(), contacts.d2,dim*sizeof(*contacts.d2));
    return make_tuple(iarr, jarr, darr);
}

namespace {
    // exclusions of (Id a, Id b) pairs; include both (a, b) and (b, a) in
    // the high and low order bits.
    typedef std::unordered_set<uint64_t> SpatialHashExclusions;

    struct Exclusions : SpatialHashExclusions {
        void add(Id i, Id j) {
            uint64_t key1(i), key2(j);
            key1 <<= 32;
            key2 <<= 32;
            key1 |= j;
            key2 |= i;
            insert(key1);
            insert(key2);
        }
    };
}

static object hash_find_pairlist(SpatialHash& hash,
                                 float r,
                                 Exclusions const& excl,
                                 bool reuse_voxels) {

    SpatialHash::contact_array_t contacts;
    if (!reuse_voxels) {
        hash.voxelize(r);
    }
    hash.findPairlistReuseVoxels(r, excl, &contacts);

    ssize_t dim = contacts.count;
    std::for_each(contacts.d2, contacts.d2+dim, [](float& x) {x=std::sqrt(x);});

    auto iarr = ids_t(dim);
    auto jarr = ids_t(dim);
    auto darr = pos_t(dim);
    memcpy(iarr.mutable_data(), contacts.i, dim*sizeof(*contacts.i));
    memcpy(jarr.mutable_data(), contacts.j, dim*sizeof(*contacts.j));
    memcpy(darr.mutable_data(), contacts.d2,dim*sizeof(*contacts.d2));
    return make_tuple(iarr, jarr, darr);
}

namespace desres { namespace msys {
    void export_spatial_hash(module m) {

        class_<Exclusions>(m, "SpatialHashExclusions")
            .def(init<>())
            .def("add", &Exclusions::add)
            ;

        class_<SpatialHash>(m, "SpatialHash")
            .def(init([](pos_t posarr, ids_t idsarr, object boxobj) {
                if (posarr.ndim() != 2 || posarr.shape(1) != 3) {
                    throw std::invalid_argument("expected Nx3 array for pos");
                }
                auto boxarr = box_t::ensure(boxobj);
                auto box = boxobj.is_none() ? nullptr : boxarr.data();
                if (box && (boxarr.ndim() != 2 || boxarr.shape(0)!=3 || boxarr.shape(1) != 3)) {
                    throw std::invalid_argument("expected 3x3 array or none for box");
                }
                return new SpatialHash(posarr.data(), idsarr.shape(0), idsarr.data(), box);
                }), arg("pos"), arg("ids"), arg("box")=none())
            .def("voxelize", &SpatialHash::voxelize, return_value_policy::reference)
            .def("findWithin", hash_find_within,
                     arg("r"), 
                     arg("pos"), 
                     arg("ids"), 
                     arg("reuse_voxels")=false)
            .def("findNearest", hash_find_nearest,
                     arg("k"), 
                     arg("pos"), 
                     arg("ids"))
            .def("findContacts", hash_find_contacts,
                     arg("r"),
                     arg("pos"),
                     arg("ids")=none(),
                     arg("reuse_voxels")=false)
            .def("findPairlist", hash_find_pairlist,
                     arg("r"),
                     arg("excl"),
                     arg("reuse_voxels")=false)
            ;
    }
}}

