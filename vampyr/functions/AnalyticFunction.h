#pragma once

#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <MRCPP/functions/AnalyticFunction.h>
//#include <MRCPP/functions/RepresentableFunction.h>

#include "PyRepresentableFunction.h"

namespace vampyr {
template <int D> void analytic_function(pybind11::module &m) {
    using namespace mrcpp;
    namespace py = pybind11;
    using namespace pybind11::literals;

    py::class_<AnalyticFunction<D>, RepresentableFunction<D>>(m, "AnalyticFunction")
        .def(py::init<std::function<double(const Coord<D> &r)>,
                      const std::vector<double> &,
                      const std::vector<double> &>(),
             "f"_a,
             "a"_a = std::vector<double>(),
             "b"_a = std::vector<double>());
}
} // namespace vampyr
