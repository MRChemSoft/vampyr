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

    py::class_<DerivativeOperator<D>>(m,
                                      "DerivativeOperator",
                                      // clang-format off
    R"mydelimiter(
        An abstract base class for derivative operators
    )mydelimiter")
        // clang-format on
        .def(py::init<const MultiResolutionAnalysis<D> &, int, int>(), "mra"_a, "root"_a, "reach"_a)
        .def("getOrder", &DerivativeOperator<D>::getOrder)
        .def(
            "__call__",
            [](DerivativeOperator<D> &oper, FunctionTree<D, double> *inp, int axis) {
                auto out = std::make_unique<FunctionTree<D, double>>(inp->getMRA());
                apply(*out, oper, *inp, axis);
                return out;
            },
            "inp"_a,
            "axis"_a = 0);

    py::class_<ABGVOperator<D>, DerivativeOperator<D>>(m,
                                                       "ABGVDerivative",
                                                       // clang-format off
    R"mydelimiter(
        ABGV derivative operator. :cite:`Alpert2002-sr`.
    )mydelimiter")
        // clang-format on
        .def(py::init<const MultiResolutionAnalysis<D> &, double, double>(), "mra"_a, "a"_a, "b"_a);

    py::class_<PHOperator<D>, DerivativeOperator<D>>(m, "PHDerivative")
        .def(py::init<const MultiResolutionAnalysis<D> &, int>(), "mra"_a, "order"_a = 1);

    py::class_<BSOperator<D>, DerivativeOperator<D>>(m,
                                                     "BSDerivative",
                                                     // clang-format off
    R"mydelimiter(
        B-spline derivative operator. :cite:`Anderson2019-bx`.
    )mydelimiter")
        // clang-format on
        .def(py::init<const MultiResolutionAnalysis<D> &, int>(), "mra"_a, "order"_a = 1);


    // Factory function to create derivative operators based on type
    m.def("Derivative", [](const MultiResolutionAnalysis<D> &mra, const std::string &type = "center", int order = 1) -> std::unique_ptr<DerivativeOperator<D>> {
        if (type == "center") {
            return std::make_unique<ABGVOperator<D>>(mra, 0.5, 0.5);
        } else if (type == "simple") {
            return std::make_unique<ABGVOperator<D>>(mra, 0.0, 0.0);
        } else if (type == "forward") {
            return std::make_unique<ABGVOperator<D>>(mra, 0.0, 1.0);
        } else if (type == "backward") {
            return std::make_unique<ABGVOperator<D>>(mra, 1.0, 0.0);
        } else if (type == "b-spline") {
            return std::make_unique<BSOperator<D>>(mra, order);
        } else if (type == "ph") {
            return std::make_unique<PHOperator<D>>(mra, order);
        } else {
            throw std::invalid_argument("Unknown derivative type: " + type);
        }
    }, "mra"_a, "type"_a, "order"_a = 1);

}
} // namespace vampyr
