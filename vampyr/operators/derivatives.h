#pragma once

#include <pybind11/pybind11.h>

#include <MRCPP/operators/ABGVOperator.h>
#include <MRCPP/operators/PHOperator.h>
#include <MRCPP/operators/BSOperator.h>

namespace vampyr {

template <int D> void derivatives(pybind11::module &m) {
    using namespace mrcpp;
    namespace py = pybind11;
    using namespace pybind11::literals;

    py::class_<DerivativeOperator<D>>(m, "DerivativeOperator")
        .def(py::init<const MultiResolutionAnalysis<D> &>(), "mra"_a)
        .def("getOrder", &ABGVOperator<D>::getOrder);

    py::class_<ABGVOperator<D>, DerivativeOperator<D>>(m, "ABGVOperator")
        .def(py::init<const MultiResolutionAnalysis<D> &, double, double>(),
             "mra"_a,
             "a"_a,
             "b"_a);

    py::class_<PHOperator<D>, DerivativeOperator<D>>(m, "PHOperator")
        .def(py::init<const MultiResolutionAnalysis<D> &, int>(),
             "mra"_a,
             "order"_a = 1);

    py::class_<BSOperator<D>, DerivativeOperator<D>>(m, "BSOperator")
        .def(py::init<const MultiResolutionAnalysis<D> &, int>(),
             "mra"_a,
             "order"_a = 1);

}
} // namespace vampyr
