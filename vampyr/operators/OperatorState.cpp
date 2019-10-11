#include "pybind11/pybind11.h"

#include "trees/MWNode.h"

#include "operators/OperatorState.h"

using namespace mrcpp;
namespace py = pybind11;

namespace vampyr {

void operator_state(py::module &m) {
    const auto D = 3;

// The double * wont work properly in vampyr since raw pointers does not exist in python. 
// This should probably be changed to an std::vector in mrcpp. 
    //OperatorState
    py::class_<OperatorState<D>> operatorstate(m, "OperatorState");
    operatorstate.def(py::init<MWNode<D> &, double*>());


}
} // namespace vampyr
