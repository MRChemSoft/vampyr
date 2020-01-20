#pragma once

#include <pybind11/pybind11.h>

#include <MRCPP/operators/OperatorStatistics.h>
#include <MRCPP/trees/MWNode.h>
#include <MRCPP/trees/MultiResolutionAnalysis.h>

namespace vampyr {
template <int D> void operator_statistics(pybind11::module &m) {
    using namespace mrcpp;
    namespace py = pybind11;

    //OperatorStatistics
    py::class_<OperatorStatistics<D>> operatorstatistics(m, "OperatorStatistics");
    operatorstatistics.def(py::init<>())
        .def("flushNodeCounters", &OperatorStatistics<D>::flushNodeCounters, "flushNodeCounters")
        .def("incrementFNodeCounters", &OperatorStatistics<D>::incrementFNodeCounters, "incrementFNodeCounters")
        .def("incrementGNodeCounters", &OperatorStatistics<D>::incrementGNodeCounters, "incrementGNodeCounters");
}
} // namespace vampyr
