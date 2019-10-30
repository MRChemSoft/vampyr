#include <array>

#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/functional.h>
#include <pybind11/stl.h>

#include "trees/MultiResolutionAnalysis.h"

#include "operators/MWOperator.h"

using namespace mrcpp;
namespace py = pybind11;

namespace vampyr {

void mw_operator(py::module &m) {

    //MWOperator
    py::class_<MWOperator> mwoperator(m, "MWOperator");
    mwoperator.def(py::init<MultiResolutionAnalysis<2>>())
        .def("size", &MWOperator::size, "Size")
        .def("push_back", &MWOperator::push_back, "Push")
        .def("clear", &MWOperator::clear, "Clear")
        .def("getMaxBandWidth", &MWOperator::getMaxBandWidth, "Get")
        .def("calcBandWidths", &MWOperator::calcBandWidths, "Calc")
        .def("clearBandWidths", &MWOperator::clearBandWidths, "Clear");


}
} // namespace vampyr
