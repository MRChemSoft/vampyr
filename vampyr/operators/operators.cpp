#include <array>

#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "trees/FunctionTree.h"
#include "trees/FunctionTreeVector.h"

#include "trees/MultiResolutionAnalysis.h"

#include "PyRepresentableFunction.h"

#include "operators/ConvolutionOperator.h"
#include "operators/DerivativeOperator.h"
#include "operators/GreensKernel.h"

#include "treebuilders/apply.h"
#include "treebuilders/grid.h"

#include "functions/RepresentableFunction.h"

using namespace mrcpp;
namespace py = pybind11;
using namespace pybind11::literals;

namespace vampyr {

void abgv_operator(py::module &, py::class_<DerivativeOperator<3>> &derivativeoperator);
void bs_operator(py::module &, py::class_<DerivativeOperator<3>> &derivativeoperator);
void derivative_kernel(py::module &);
void helmholtz_kernel(py::module &, py::class_<GreensKernel> &greenskernel);
void helmholtz_operator(py::module &);
void identityconvolution_operator(py::module &);
void identity_kernel(py::module &, py::class_<GreensKernel> &greenskernel);
void mw_operator(py::module &);
void operator_state(py::module &);
void operator_statistics(py::module &);
void ph_operator(py::module &, py::class_<DerivativeOperator<3>> &derivativeoperator);
void poisson_operator(py::module &);

void operators(py::module &m) {

    const auto D = 3;

    //ConvolutionOperator
    py::class_<ConvolutionOperator<D>> convop(m, "ConvolutionOperator");
    convop.def(py::init<MultiResolutionAnalysis<D> &, double>());

    //DerivativeOperator
    py::class_<DerivativeOperator<D>> derivativeoperator(m, "Derivative Operator");
    derivativeoperator.def(py::init<MultiResolutionAnalysis<D>>());

    //GreensKernel
    py::class_<GreensKernel> greenskernel(m, "GreensKernel");

    m.def(
        "divergence",
        py::overload_cast<FunctionTree<D> &, DerivativeOperator<D> &, FunctionTreeVector<D> &>(&mrcpp::divergence<D>));

    m.def("gradient", py::overload_cast<DerivativeOperator<D> &, FunctionTree<D> &>(&gradient<D>));

    abgv_operator(m, derivativeoperator);
    bs_operator(m, derivativeoperator);
    derivative_kernel(m);
    helmholtz_kernel(m, greenskernel);
    helmholtz_operator(m);
    identityconvolution_operator(m);
    identity_kernel(m, greenskernel);
    mw_operator(m);
    operator_state(m);
    operator_statistics(m);
    ph_operator(m, derivativeoperator);
    poisson_operator(m);
}
} // namespace vampyr
