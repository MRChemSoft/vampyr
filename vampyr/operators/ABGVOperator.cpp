#include <array>

#include "pybind11/eigen.h"
#include "pybind11/pybind11.h"
#include "pybind11/functional.h"
#include "pybind11/stl.h"

#include "trees/MultiResolutionAnalysis.h"

#include "operators/DerivativeOperator.h"
#include "operators/ABGVOperator.h"

using namespace mrcpp;
namespace py = pybind11;

namespace vampyr {

void abgv_operator(py::module &m, py::class_<DerivativeOperator<3>> &derivativeoperator) {
    const auto D = 3;


   //ABGVOperator
    py::class_<ABGVOperator<D>> abgvoperator(m, "ABGVOperator", derivativeoperator);
    abgvoperator.def(py::init<MultiResolutionAnalysis<D> &, double, double>(),
             "Derivative Operator: The ABGVOperator for differentiation");


}
} // namespace vampyr
