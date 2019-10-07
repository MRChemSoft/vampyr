#include "pybind11/pybind11.h"

#include "trees/MultiResolutionAnalysis.h"

#include "operators/ConvolutionOperator.h"

#include "functions/RepresentableFunction.h"

using namespace mrcpp;
namespace py = pybind11;

namespace vampyr {

void convolution_operator(py::module &m) {
    const auto D = 3;

    //ConvolutionOperator
	  py::class_<ConvolutionOperator<D>> convop(m, "ConvolutionOperator");
      convop.def(py::init<MultiResolutionAnalysis<D> &, double>());


}
} // namespace vampyr
