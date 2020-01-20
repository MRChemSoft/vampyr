#pragma once

#include <pybind11/pybind11.h>

#include <MRCPP/operators/DerivativeOperator.h>
#include <MRCPP/operators/PHOperator.h>
#include <MRCPP/trees/MultiResolutionAnalysis.h>

namespace vampyr {
template <int D>
void ph_operator(pybind11::module &m, pybind11::class_<mrcpp::DerivativeOperator<D>> &derivativeoperator) {
    using namespace mrcpp;
    namespace py = pybind11;

    //PHOperator
    py::class_<PHOperator<D>>(m, "PHOperator", derivativeoperator)
        .def(py::init<const MultiResolutionAnalysis<D> &, int>());
}
} // namespace vampyr
