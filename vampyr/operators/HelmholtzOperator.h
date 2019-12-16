#pragma once

#include <pybind11/pybind11.h>

#include <MRCPP/operators/ConvolutionOperator.h>
#include <MRCPP/operators/HelmholtzOperator.h>
#include <MRCPP/trees/MultiResolutionAnalysis.h>

namespace vampyr {
template <int D> void helmholtz_operator(pybind11::module &) {}

template <> void helmholtz_operator<1>(pybind11::module &) {}
template <> void helmholtz_operator<2>(pybind11::module &) {}
template <> void helmholtz_operator<3>(pybind11::module &m) {
    using namespace mrcpp;
    namespace py = pybind11;
    using namespace pybind11::literals;

    //HelmholtzOperator
    py::class_<HelmholtzOperator, ConvolutionOperator<3>>(m, "HelmholtzOperator")
        .def(py::init<const MultiResolutionAnalysis<3> &, double, double>(),
             "MRA"_a,
             "mu"_a,
             "precision"_a = -1.0,
             "ConvolutionOperator: HelmholtzOperator exp(-mu*r)/|r-r'|");
}
} // namespace vampyr
