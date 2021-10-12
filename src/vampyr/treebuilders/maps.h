#pragma once

#include <pybind11/functional.h>

#include "PyFunctionMap.h"

namespace vampyr {
template <int D> void map(pybind11::module &m) {
    using namespace mrcpp;
    namespace py = pybind11;
    using namespace pybind11::literals;

    py::class_<PyFunctionMap<D>>(m, "FunctionMap")
        .def(py::init<std::function<double (double)>, double>(),
            "fmap"_a,
            "prec"_a)
        .def("__call__", [](PyFunctionMap<D> &F, FunctionTree<D> &inp){
                auto old_threads = mrcpp_get_num_threads();
                set_max_threads(1);
                auto out = F(inp);
                set_max_threads(old_threads);
                return out;
            },
            "inp"_a);
}

template <int D> void advanced_map(pybind11::module &m) {
    using namespace mrcpp;
    namespace py = pybind11;
    using namespace pybind11::literals;

    m.def("map", [] (double prec,
                     FunctionTree<D> &out,
                     FunctionTree<D> &inp,
                     std::function<double(double)> fmap,
                     int max_iter,
                     bool abs_prec) {
              auto old_threads = mrcpp_get_num_threads();
              mrcpp::set_max_threads(1);
              mrcpp::map<D>(prec, out, inp, fmap, max_iter, abs_prec);
              mrcpp::set_max_threads(old_threads);
          },
          "prec"_a = -1.0,
          "out"_a,
          "inp"_a,
          "fmap"_a,
          "max_iter"_a = -1,
          "abs_prec"_a = false);
}
} // namespace vampyr
