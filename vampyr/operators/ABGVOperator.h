#pragma once

#include <pybind11/pybind11.h>

#include <MRCPP/operators/ABGVOperator.h>
#include <MRCPP/operators/DerivativeOperator.h>
#include <MRCPP/trees/MultiResolutionAnalysis.h>

namespace vampyr {
template <int D>
void abgv_operator(pybind11::module &m, pybind11::class_<mrcpp::DerivativeOperator<D>> &derivativeoperator) {
    using namespace mrcpp;
    namespace py = pybind11;

    //ABGVOperator
    py::class_<ABGVOperator<D>>(m, "ABGVOperator", derivativeoperator)
        .def(py::init<const MultiResolutionAnalysis<D> &, double, double>(),
             "Derivative Operator: The ABGVOperator for differentiation");
}
} // namespace vampyr
