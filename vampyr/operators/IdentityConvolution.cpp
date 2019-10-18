#include <array>

#include "pybind11/eigen.h"
#include "pybind11/pybind11.h"
#include "pybind11/functional.h"
#include "pybind11/stl.h"

#include "trees/MultiResolutionAnalysis.h"

#include "operators/ConvolutionOperator.h"
#include "operators/IdentityConvolution.h"

using namespace mrcpp;
namespace py = pybind11;

namespace vampyr {

void identityconvolution_operator(py::module &m, py::class_<ConvolutionOperator<3>> &convop) {
    const auto D = 3;

    
    //IdentityConvolution
	  py::class_<IdentityConvolution<D>>identityconvop(m, "IdentityConvolution", convop);
      identityconvop.def(py::init<MultiResolutionAnalysis<D> &, double>());


}
} // namespace vampyr
