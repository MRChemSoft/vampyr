#pragma once

#include <pybind11/pybind11.h>

#include <MRCPP/operators/BSOperator.h>
#include <MRCPP/operators/DerivativeOperator.h>
#include <MRCPP/trees/MultiResolutionAnalysis.h>

namespace vampyr {
template <int D>
void bs_operator(pybind11::module &m, pybind11::class_<mrcpp::DerivativeOperator<D>> &derivativeoperator) {
    using namespace mrcpp;
    namespace py = pybind11;

    //BSOperator
    py::class_<BSOperator<D>>(m, "BSOperator", derivativeoperator)
        .def(py::init<const MultiResolutionAnalysis<D> &, int>());
}
} // namespace vampyr
