#include "pybind11/pybind11.h"

#include "trees/MultiResolutionAnalysis.h"
#include "trees/MWNode.h"

#include "operators/OperatorStatistics.h"

using namespace mrcpp;
namespace py = pybind11;

namespace vampyr {

void operator_statistics(py::module &m) {
    const auto D = 3;


    //OperatorStatistics
    py::class_<OperatorStatistics<D>> operatorstatistics(m, "OperatorStatistics");
    operatorstatistics.def(py::init<>())
        .def("flushNodeCounters", &OperatorStatistics<D>::flushNodeCounters, "flushNodeCounters")
        .def("incrementFNodeCounters", &OperatorStatistics<D>::incrementFNodeCounters, "incrementFNodeCounters")
        .def("incrementGNodeCounters", &OperatorStatistics<D>::incrementGNodeCounters, "incrementGNodeCounters");


}
} // namespace vampyr