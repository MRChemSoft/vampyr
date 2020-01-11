#pragma once

#include <pybind11/pybind11.h>

#include <MRCPP/functions/RepresentableFunction.h>
#include <MRCPP/operators/ConvolutionOperator.h>
#include <MRCPP/operators/DerivativeOperator.h>
#include <MRCPP/treebuilders/apply.h>
#include <MRCPP/treebuilders/grid.h>
#include <MRCPP/trees/FunctionTree.h>
#include <MRCPP/trees/FunctionTreeVector.h>
#include <MRCPP/trees/MultiResolutionAnalysis.h>

#include "ABGVOperator.h"
#include "BSOperator.h"
#include "HelmholtzOperator.h"
#include "IdentityConvolution.h"
#include "OperatorState.h"
#include "OperatorStatistics.h"
#include "PHOperator.h"
#include "PoissonOperator.h"

namespace vampyr {
void helmholtz_kernel(pybind11::module &);
void derivative_kernel(pybind11::module &);
void identity_kernel(pybind11::module &);

template <int D> void operators(pybind11::module &m) {
    using namespace mrcpp;
    namespace py = pybind11;

    //ConvolutionOperator
    py::class_<ConvolutionOperator<D>>(m, "ConvolutionOperator")
        .def(py::init<const MultiResolutionAnalysis<D> &, double>());

    identityconvolution_operator<D>(m);
    helmholtz_operator<D>(m);
    poisson_operator<D>(m);

    //DerivativeOperator
    py::class_<DerivativeOperator<D>> derivativeoperator(m, "Derivative Operator");
    derivativeoperator.def(py::init<const MultiResolutionAnalysis<D> &>());

    abgv_operator<D>(m, derivativeoperator);
    bs_operator<D>(m, derivativeoperator);
    ph_operator<D>(m, derivativeoperator);

    operator_state<D>(m);
    operator_statistics<D>(m);
}
} // namespace vampyr
