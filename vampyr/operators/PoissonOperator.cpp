#include "pybind11/pybind11.h"

#include "trees/MultiResolutionAnalysis.h"

#include "PyRepresentableFunction.h"

#include "operators/ConvolutionOperator.h"
#include "operators/PoissonOperator.h"
#include "operators/MWOperator.h"

using namespace mrcpp;
namespace py = pybind11;
using namespace pybind11::literals;

namespace vampyr {

void poisson_operator(py::module &m) {
    const auto D = 3;

    //PoissonOperator
    py::class_<PoissonOperator, ConvolutionOperator<D>, MWOperator>(m, "PoissonOperator", py::multiple_inheritance())
        .def(py::init<MultiResolutionAnalysis<D> &, double>(),
             "MRA"_a,
             "precision"_a,
             "ConvolutionOperator: PoissonOperator 1/|r-r'|");


}
} // namespace vampyr