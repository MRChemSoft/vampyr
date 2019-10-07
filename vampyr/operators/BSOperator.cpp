#include "pybind11/pybind11.h"

#include "trees/MultiResolutionAnalysis.h"

#include "operators/BSOperator.h"

using namespace mrcpp;
namespace py = pybind11;

namespace vampyr {

void bs_operator(py::module &m) {
    const auto D = 3;

    //BSOperator
 	  py::class_<BSOperator<D>> bsoperator(m, "BSOperator");
      bsoperator.def(py::init<MultiResolutionAnalysis<D> &, int>());


}
} // namespace vampyr
