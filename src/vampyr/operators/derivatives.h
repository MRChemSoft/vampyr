#pragma once

#include <pybind11/pybind11.h>

#include <MRCPP/operators/ABGVOperator.h>
#include <MRCPP/operators/BSOperator.h>
#include <MRCPP/operators/PHOperator.h>

namespace vampyr {

template <int D> void derivatives(pybind11::module &m) {
    using namespace mrcpp;
    namespace py = pybind11;
    using namespace pybind11::literals;

    py::class_<DerivativeOperator<D>>(m, "DerivativeOperator")
        .def(py::init<const MultiResolutionAnalysis<D> &>(), "mra"_a)
        .def("getOrder", &DerivativeOperator<D>::getOrder)
        .def(
            "__call__",
            [](DerivativeOperator<D> &oper, FunctionTree<D> *inp, int axis) {
                auto out = std::make_unique<FunctionTree<D>>(inp->getMRA());
                apply(*out, oper, *inp, axis);
                return out;
            },
            "inp"_a,
            "axis"_a = 0);

    py::class_<ABGVOperator<D>, DerivativeOperator<D>>(m, "ABGVDerivative")
        .def(py::init<const MultiResolutionAnalysis<D> &, double, double>(), "mra"_a, "a"_a, "b"_a);

    py::class_<PHOperator<D>, DerivativeOperator<D>>(m, "PHDerivative")
        .def(py::init<const MultiResolutionAnalysis<D> &, int>(), "mra"_a, "order"_a = 1);

    py::class_<BSOperator<D>, DerivativeOperator<D>>(m, "BSDerivative")
        .def(py::init<const MultiResolutionAnalysis<D> &, int>(), "mra"_a, "order"_a = 1);
}
} // namespace vampyr
