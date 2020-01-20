#pragma once

#include <pybind11/pybind11.h>

#include <MRCPP/operators/ConvolutionOperator.h>
#include <MRCPP/operators/IdentityConvolution.h>
#include <MRCPP/trees/MultiResolutionAnalysis.h>

namespace vampyr {
template <int D> void identityconvolution_operator(pybind11::module &m) {
    using namespace mrcpp;
    namespace py = pybind11;

    //IdentityConvolution
    py::class_<IdentityConvolution<D>, ConvolutionOperator<D>>(m, "IdentityConvolution")
        .def(py::init<const MultiResolutionAnalysis<D> &, double>());
}
} // namespace vampyr
