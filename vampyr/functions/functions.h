#pragma once

#include "PyRepresentableFunction.h"

namespace vampyr {

template <int D> void functions(pybind11::module &m) {
    using namespace mrcpp;
    namespace py = pybind11;
    using namespace pybind11::literals;

    // RepresentableFunction class
    py::class_<RepresentableFunction<D>, PyRepresentableFunction<D>>(m, "RepresentableFunction")
        .def(py::init<const Coord<D> &, const Coord<D> &>())
        .def("__call__", [](const RepresentableFunction<D> &func, const Coord<D> &r) {
            return func.evalf(r);
        });

}
} // namespace vampyr
