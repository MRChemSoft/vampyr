#pragma once

#include "PyRepresentableFunction.h"

#include "gaussians.h"

namespace vampyr {

template <int D> void functions(pybind11::module &m) {
    using namespace mrcpp;
    namespace py = pybind11;
    using namespace pybind11::literals;

    // RepresentableFunction class
    py::class_<RepresentableFunction<D, double>, PyRepresentableFunction<D>>(m,
                                                                     "RepresentableFunction",
                                                                     R"mydelimiter(
       Filler text
    )mydelimiter")
        .def(py::init<const std::vector<double> &, const std::vector<double> &>())
        .def("__call__", [](const RepresentableFunction<D, double> &func, const Coord<D> &r) { return func.evalf(r); });

    gaussians<D>(m);
}
} // namespace vampyr
