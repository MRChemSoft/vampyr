#include "pybind11/pybind11.h"

#include "trees/MultiResolutionAnalysis.h"

#include "operators/ABGVOperator.h"

using namespace mrcpp;
namespace py = pybind11;

namespace vampyr {

void abgv_operator(py::module &m) {
    const auto D = 3;


   //ABGVOperator
    py::class_<ABGVOperator<D>> abgvoperator(m, "ABGVOperator");
    abgvoperator.def(py::init<MultiResolutionAnalysis<D> &, double, double>(),
             "Derivative Operator: The ABGVOperator for differentiation");


}
} // namespace vampyr