#pragma once

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <MRCPP/operators/OperatorState.h>
#include <MRCPP/trees/MWNode.h>

namespace vampyr {
template <int D> void operator_state(pybind11::module &m) {
    using namespace mrcpp;
    namespace py = pybind11;

    py::class_<OperatorState<D>>(m, "OperatorState").def(py::init<MWNode<D> &, std::vector<double>>());
}
} // namespace vampyr
