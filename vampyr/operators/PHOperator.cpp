#include <array>

#include "pybind11/eigen.h"
#include "pybind11/pybind11.h"
#include "pybind11/functional.h"
#include "pybind11/stl.h"

#include "trees/MultiResolutionAnalysis.h"

#include "operators/PHOperator.h"

using namespace mrcpp;
namespace py = pybind11;

namespace vampyr {

void ph_operator(py::module &m) {
    const auto D = 3;

    //PHOperator
    py::class_<PHOperator<D>> phoperator(m, "PHOperator");
    phoperator.def(py::init<MultiResolutionAnalysis<D> &, int>());


}
} // namespace vampyr
