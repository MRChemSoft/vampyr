#include "pybind11/pybind11.h"

#include "trees/MultiResolutionAnalysis.h"

#include "operators/DerivativeOperator.h"

using namespace mrcpp;
namespace py = pybind11;

namespace vampyr {

void derivative_operator(py::module &m) {
    const auto D = 3;

    //DerivativeOperator
    py::class_<DerivativeOperator<D>> derivativeoperator(m, "Derivative Operator");
    derivativeoperator.def(py::init<MultiResolutionAnalysis<D>>());


}
} // namespace vampyr
