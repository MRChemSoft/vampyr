#include "pybind11/pybind11.h"

#include "trees/MultiResolutionAnalysis.h"

#include "operators/ConvolutionOperator.h"
#include "operators/HelmholtzOperator.h"

#include "treebuilders/apply.h"

using namespace mrcpp;
namespace py = pybind11;
using namespace pybind11::literals;

namespace vampyr {

void helmholtz_operator(py::module &m) {
    const auto D = 3;

    //HelmholtzOperator
    py::class_<HelmholtzOperator, ConvolutionOperator<D>>(m, "HelmholtzOperator", py::multiple_inheritance())
        .def(py::init<MultiResolutionAnalysis<D> &, double, double>(),
             "MRA"_a,
             "mu"_a,
             "precision"_a = -1.0,
             "ConvolutionOperator: HelmholtzOperator exp(-mu*r)/|r-r'|");


}
} // namespace vampyr
