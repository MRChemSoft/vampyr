#pragma once

#include <pybind11/pybind11.h>

#include <MRCPP/operators/OperatorState.h>
#include <MRCPP/trees/MWNode.h>

namespace vampyr {
template <int D> void operator_state(pybind11::module &m) {
    using namespace mrcpp;
    namespace py = pybind11;

    // The double * wont work properly in vampyr since raw pointers does not exist in python.
    // This should probably be changed to an std::vector in mrcpp.
    //OperatorState
    py::class_<OperatorState<D>>(m, "OperatorState").def(py::init<MWNode<D> &, double *>());
}
} // namespace vampyr
