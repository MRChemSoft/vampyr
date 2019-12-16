#pragma once

#include <pybind11/pybind11.h>

#include <MRCPP/operators/ConvolutionOperator.h>
#include <MRCPP/operators/PoissonOperator.h>
#include <MRCPP/trees/MultiResolutionAnalysis.h>

namespace vampyr {
template <int D> void poisson_operator(pybind11::module &) {}

template <> void poisson_operator<1>(pybind11::module &) {}
template <> void poisson_operator<2>(pybind11::module &) {}
template <> void poisson_operator<3>(pybind11::module &m) {
    namespace py = pybind11;
    using namespace mrcpp;
    using namespace pybind11::literals;

    //PoissonOperator
    py::class_<PoissonOperator, ConvolutionOperator<3>>(m, "PoissonOperator")
        .def(py::init<const MultiResolutionAnalysis<3> &, double>(),
             "MRA"_a,
             "precision"_a,
             "ConvolutionOperator: PoissonOperator 1/|r-r'|");
}
} // namespace vampyr
